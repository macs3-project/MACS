/* Test-only allocator for an isolated build of PeakModel.py. */
#include <pthread.h>
#include <stddef.h>
#include <stdlib.h>

static pthread_mutex_t pm_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t pm_condition = PTHREAD_COND_INITIALIZER;
static pthread_t pm_main;
static void *pm_owned[4096];
static size_t pm_calls, pm_target, pm_failures, pm_live, pm_errors;
static size_t pm_failed_size, pm_failed_growth, pm_after_failure;
static int pm_worker_only, pm_failed_in_worker;

static void pm_reset(size_t target, int worker_only) {
    pthread_mutex_lock(&pm_lock);
    /* Never erase evidence of a leak or invalid free between runs. */
    pm_calls = pm_failures = pm_failed_size = pm_failed_growth = 0;
    pm_after_failure = 0;
    pm_failed_in_worker = 0;
    pm_target = target;
    pm_worker_only = worker_only;
    pm_main = pthread_self();
    pthread_mutex_unlock(&pm_lock);
}

static void *pm_realloc(void *ptr, size_t size) {
    size_t slot;
    void *result;
    pthread_mutex_lock(&pm_lock);
    /* Hold the caller's scan until a worker reaches the injected failure.
       The two-task test consequently checks worker exception propagation
       even on machines where the caller wins the scheduling race. */
    while (pm_worker_only && !pm_failures &&
           pthread_equal(pthread_self(), pm_main))
        pthread_cond_wait(&pm_condition, &pm_lock);
    pm_calls++;
    if (pm_failures)
        pm_after_failure++;
    if (pm_target && pm_calls == pm_target) {
        pm_failures++;
        pm_failed_size = size;
        pm_failed_growth = ptr != NULL;
        pm_failed_in_worker = !pthread_equal(pthread_self(), pm_main);
        pthread_cond_broadcast(&pm_condition);
        pthread_mutex_unlock(&pm_lock);
        return NULL;
    }
    for (slot = 0; slot < 4096; slot++)
        if (ptr ? pm_owned[slot] == ptr : pm_owned[slot] == NULL)
            break;
    if (slot == 4096 || size == 0) {
        pm_errors++;
        pthread_mutex_unlock(&pm_lock);
        return NULL;
    }
    result = realloc(ptr, size);
    if (result) {
        pm_owned[slot] = result;
        if (!ptr)
            pm_live++;
    }
    pthread_mutex_unlock(&pm_lock);
    return result;
}

static void pm_free(void *ptr) {
    size_t slot;
    if (!ptr)
        return;
    pthread_mutex_lock(&pm_lock);
    for (slot = 0; slot < 4096; slot++)
        if (pm_owned[slot] == ptr)
            break;
    if (slot == 4096) {
        pm_errors++;
    } else {
        free(ptr);
        pm_owned[slot] = NULL;
        pm_live--;
    }
    pthread_mutex_unlock(&pm_lock);
}

static void pm_stats(size_t *out) {
    pthread_mutex_lock(&pm_lock);
    out[0] = pm_calls;
    out[1] = pm_failures;
    out[2] = pm_live;
    out[3] = pm_errors;
    out[4] = pm_failed_size;
    out[5] = pm_failed_growth;
    out[6] = pm_failed_in_worker;
    out[7] = pm_after_failure;
    pthread_mutex_unlock(&pm_lock);
}
