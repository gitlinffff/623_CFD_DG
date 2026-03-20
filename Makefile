BUILD := build

.PHONY: all run clean test

# Configure cmake once (only if build/Makefile doesn't exist yet).
# Incremental rebuilds just call make inside build/.
$(BUILD)/Makefile:
	mkdir -p $(BUILD)
	cd $(BUILD) && cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=OFF

all: $(BUILD)/Makefile
	@$(MAKE) -C $(BUILD) main --no-print-directory

run: all
	@$(BUILD)/main

test: $(BUILD)/Makefile
	@cd $(BUILD) && cmake .. -DBUILD_TESTS=ON > /dev/null 2>&1 || true
	@$(MAKE) -C $(BUILD) test_freestream test_quad --no-print-directory
	@echo "--- test_quad ---"
	@$(BUILD)/test/test_quad
	@echo "--- test_freestream ---"
	@$(BUILD)/test/test_freestream

clean:
	rm -rf $(BUILD) results
