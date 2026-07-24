import type {
  AnalysisData,
  AnalysisRequest,
  DemoResponse,
  SequenceRecord,
} from '../types'

async function responseError(response: Response): Promise<Error> {
  try {
    const body = (await response.json()) as { detail?: string }
    return new Error(body.detail ?? `Request failed with status ${response.status}.`)
  } catch {
    return new Error(`Request failed with status ${response.status}.`)
  }
}

export async function loadDemo(): Promise<DemoResponse> {
  const response = await fetch('/api/demo')
  if (!response.ok) throw await responseError(response)
  return (await response.json()) as DemoResponse
}

export async function analyze(request: AnalysisRequest): Promise<AnalysisData> {
  const response = await fetch('/api/analyze', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(request),
  })
  if (!response.ok) throw await responseError(response)
  return (await response.json()) as AnalysisData
}

export async function parseFasta(file: File): Promise<SequenceRecord> {
  const form = new FormData()
  form.append('file', file)
  const response = await fetch('/api/records/parse', {
    method: 'POST',
    body: form,
  })
  if (!response.ok) throw await responseError(response)
  return (await response.json()) as SequenceRecord
}

export async function downloadExport(
  kind: 'json' | 'csv',
  request: AnalysisRequest,
): Promise<void> {
  const response = await fetch(`/api/export/${kind}`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(request),
  })
  if (!response.ok) throw await responseError(response)

  const blob = await response.blob()
  const disposition = response.headers.get('Content-Disposition') ?? ''
  const match = disposition.match(/filename="?([^";]+)"?/i)
  const filename = match?.[1] ?? `genome_embeddings_analysis.${kind}`
  const url = URL.createObjectURL(blob)
  const anchor = document.createElement('a')
  anchor.href = url
  anchor.download = filename
  document.body.append(anchor)
  anchor.click()
  anchor.remove()
  URL.revokeObjectURL(url)
}
