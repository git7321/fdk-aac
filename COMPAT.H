/*
 * Copyright (C) nu774
 * For conditions of distribution and use, see copyright notice in COPYING
 */
#ifndef COMPAT_H
#define COMPAT_H

int64_t aacenc_timer(void);
FILE *aacenc_fopen(const char *name, const char *mode);
#ifdef _WIN32
void aacenc_getmainargs(int *argc, char ***argv);
#else
#define aacenc_getmainargs(argc, argv) (void)0
#endif
char *aacenc_to_utf8(const char *s);
#ifdef _WIN32
int aacenc_fprintf(FILE *fp, const char *fmt, ...);
#else
#define aacenc_fprintf fprintf
#endif
const char *aacenc_basename(const char *path);
int aacenc_seekable(FILE *fp);

#endif
