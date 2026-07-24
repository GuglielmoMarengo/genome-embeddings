import app


def test_app_module_imports_without_starting_server():
    assert callable(app.main)
