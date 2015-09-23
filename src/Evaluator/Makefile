LS_THOMAS_LEGACY = -DLS_REPEAT_FIRST_CITY_AT_END \
	-DLS_SOLUTION_LONGINT -DLS_DISTANCE_LONGINT \
	-DTESTLS_SOLUTION_LONGINT -DTESTLS_DISTANCE_LONGINT

LS_LEONORA_LEGACY = -DLS_SOLUTION_LONGINT -DLS_DISTANCE_DOUBLE \
	-DTESTLS_SOLUTION_LONGINT -DTESTLS_DISTANCE_DOUBLE


GSLFLAGS = -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF

CPPFLAGS = $(GSLFLAGS) $(LS_LEONORA_LEGACY)

# For deactivating assertiond, add  `-DNDEBUG' to the CPPFLAGS:
CPPFLAGS += -DNDEBUG 

# For extra debug output, add `-DLS_EXTRA_DEBUG_OUTPUT' to the CPPFLAGS:
#CPPFLAGS +=-DLS_EXTRA_DEBUG_OUTPUT

# For obtaining the `array' version of the library, add `-DLS_ARRAY'
# to the CPPFLAGS
CFLAGS = -std=gnu99 -Wall -Winline -pedantic
CFLAGS += -O3 
#CFLAGS += -ggdb

LDFLAGS = -lgsl -lgslcblas

# Uncomment the following two lines for including profiling information
#CFLAGS += -pg      
#LDFLAGS += -pg

LS_SOURCES = evaluate.c 

SOURCES = evaluate.c $(LS_SOURCES) 
EXECUTABLES = evaluate

all: $(EXECUTABLES)

$(EXECUTABLES): evaluate% : $(SOURCES:.c=%.o) 
	@echo -e "Linking $@ from:"; \
	 for i in $^; do echo -e "\t\t$$i"; done
	@$(CC) $(LDFLAGS) -o $@ $^ 

clean:
	@echo -e "Cleaning directory from *.o and *.d files"
	@rm -f *.o *.d gmon.out

distclean: clean
	@echo "Removing target"
	@rm -f $(EXECUTABLES)

ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
ifneq ($(findstring tarball,$(MAKECMDGOALS)),tarball)
-include $(SOURCES:%.c=%.d)
endif
endif

tarball: distclean
	@echo "Creating tarball"
	@cd ..; tar cf Evaluator.tgz Evaluator


%.d: %.c
	@printf "Making dependencies for %-20s\n" $<
	@$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.o: %.c  | Makefile
	@printf "Compiling file %-29s\n" $<
	@$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

.PHONY: all clean dist-clean tarball
