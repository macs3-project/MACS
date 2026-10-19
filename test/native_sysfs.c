/* Test-only Linux discovery shim. Redirect only the allocator's two sysfs
 * reads to fixtures; leave allocation, mmap/mremap/madvise and MACS3 intact.
 * Never changes host transparent-huge-page settings. */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int redirect_open(const char *symbol, const char *path, int flags, mode_t mode)
{
    int (*real_open)(const char *, int, ...) = dlsym(RTLD_NEXT, symbol);
    const char *root = getenv("MACS3_TEST_SYSFS");
    const char *prefix = "/sys/kernel/mm/transparent_hugepage/";
    char replacement[4096];
    if (root && path && !strncmp(path, prefix, strlen(prefix))) {
        const char *name = path + strlen(prefix);
        if (!strcmp(name, "enabled") || !strcmp(name, "hpage_pmd_size")) {
            int n = snprintf(replacement, sizeof(replacement), "%s/%s", root, name);
            if (n > 0 && (size_t)n < sizeof(replacement))
                path = replacement;
        }
    }
    return real_open(path, flags, mode);
}

static int has_mode_argument(int flags)
{
    if (flags & O_CREAT)
        return 1;
#ifdef O_TMPFILE
    if ((flags & O_TMPFILE) == O_TMPFILE)
        return 1;
#endif
    return 0;
}

#define WRAP_OPEN(name) \
int name(const char *path, int flags, ...) \
{ \
    mode_t mode = 0; \
    if (has_mode_argument(flags)) { \
        va_list args; \
        va_start(args, flags); \
        mode = va_arg(args, mode_t); \
        va_end(args); \
    } \
    return redirect_open(#name, path, flags, mode); \
}
WRAP_OPEN(open)
WRAP_OPEN(open64)
