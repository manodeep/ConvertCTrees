#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

FILE *check_fopen(char *filename, char *mode) {
  FILE *res = fopen(filename, mode);
  if (res == NULL) {
    if (mode[0] == 'w')
      fprintf(stderr, "Failed to open file %s for writing!\n", filename);
    else if (mode[0] == 'a')
      fprintf(stderr, "Failed to open file %s for appending!\n", filename);
    else
      fprintf(stderr, "Failed to open file %s for reading!\n", filename);
    exit(1);
  }
  return res;
}

void *check_realloc(void *ptr, size_t size, char *reason) {
  void *res = realloc(ptr, size);
  if ((res == NULL) && (size > 0)) {
    fprintf(stderr, "[Error] Failed to allocate %zu bytes of memory (%s)!\n", size, reason);
    assert(0);
  }
  return res;
}

/* #ifdef CHECK_FILE_IO */
/* size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream) */
/* { */
/* 	size_t nwritten; */
/* 	nwritten = fwrite(ptr, size, nmemb, stream); */
/* 	if(nwritten != nmemb) */
/* 	{ */
/* 		fprintf(stderr,"I/O error (fwrite) has occured.\n"); */
/* 		fprintf(stderr,"Instead of writing nmemb=%zu, I got nwritten = %zu ..exiting\n",nmemb,nwritten); */
/* 		exit(EXIT_FAILURE); */
/* 	} */
/* 	return nwritten; */
/* } */
/* #else */
/* size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream) */
/* { */
/* 	return fwrite(ptr, size, nmemb, stream); */

/* } */
/* #endif */
