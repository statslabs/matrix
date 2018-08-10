//
// Created by Yi Pan (Institute of Cancer and Genomic Sciences) on 10/08/2018.
//

#ifndef STATSLABS_MATRIX_ERROR_H
#define STATSLABS_MATRIX_ERROR_H

static void err_doit(int, int, const char *, va_list);

/*
 * Nonfatal error related to a system call.
 * Print a message and return.
 */
void
err_ret(const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  err_doit(1, errno, fmt, ap);
  va_end(ap);
}

/*
 * Fatal error related to a system call.
 * Print a message and terminate.
 */
void
err_sys(const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  err_doit(1, errno, fmt, ap);
  va_end(ap);
  exit(1);
}

/*
 * Nonfatal error unrelated to a system call.
 * Error code passed as explict parameter.
 * Print a message and return.
 */
void
err_cont(int error, const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  err_doit(1, error, fmt, ap);
  va_end(ap);
}

/*
 * Fatal error unrelated to a system call.
 * Error code passed as explict parameter.
 * Print a message and terminate.
 */
void
err_exit(int error, const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  err_doit(1, error, fmt, ap);
  va_end(ap);
  exit(1);
}

/*
 * Fatal error related to a system call.
 * Print a message, dump core, and terminate.
 */
void
err_dump(const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  err_doit(1, errno, fmt, ap);
  va_end(ap);
  abort();        /* dump core and terminate */
  exit(1);        /* shouldn't get here */
}

/*
 * Nonfatal error unrelated to a system call.
 * Print a message and return.
 */
void
err_msg(const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  err_doit(0, 0, fmt, ap);
  va_end(ap);
}

/*
 * Fatal error unrelated to a system call.
 * Print a message and terminate.
 */
void
err_quit(const char *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  err_doit(0, 0, fmt, ap);
  va_end(ap);
  exit(1);
}

// Print a message and return to caller
// Caller specifies "errnoflag".
static void
err_doit(int errnoflag, int error, const char *fmt, va_list ap) {

  char buf[4096];

  vsnprintf(buf, 4096 - 1, fmt, ap);
  if (errnoflag)
    snprintf(buf + strlen(buf), 4096 - strlen(buf) - 1, ": %s",
             strerror(error));
  strcat(buf, "\n");
//  if (log_to_stderr) {
  fflush(stdout);
  fputs(buf, stderr);
  fflush(stderr);
//  } else {
//    syslog(priority, "%s", buf);
//  }
}

#endif //STATSLABS_MATRIX_ERROR_H
