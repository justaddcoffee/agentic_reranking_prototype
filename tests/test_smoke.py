def test_package_importable() -> None:
    import agentic_reranking

    assert agentic_reranking.__version__ == "0.1.0"
