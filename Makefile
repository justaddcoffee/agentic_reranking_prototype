.PHONY: lint fmt typecheck test all

lint:
	uv run ruff check src tests

fmt:
	uv run ruff format src tests

typecheck:
	uv run mypy src

test:
	uv run pytest

all: fmt lint typecheck test
