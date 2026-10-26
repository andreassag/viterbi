CARGO ?= cargo

all:
	@if command -v $(CARGO) >/dev/null 2>&1; then \
		$(CARGO) build --release; \
	else \
		echo "Cargo not in PATH. Please install rustup or run scripts/test_example.sh"; \
	fi

test:
	@if command -v $(CARGO) >/dev/null 2>&1; then \
		$(CARGO) test --release; \
	else \
		python3 scripts/verify_logic.py; \
	fi

clean:
	@if command -v $(CARGO) >/dev/null 2>&1; then \
		$(CARGO) clean; \
	fi
	rm -rf target *.hmm *.tsv

.PHONY: all test clean
