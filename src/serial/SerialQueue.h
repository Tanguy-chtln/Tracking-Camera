#ifndef __QUEUE_H__
#define __QUEUE_H__

#include <stdbool.h>

#define STRING_INITIALIZER {NULL, 0}


enum buffer_queue_mode_t {
    QUEUE = 0,
    KEEP_DATA = 1
};

struct string_t {
    char *text;
    unsigned int textLength;
};

void string_init(struct string_t *string);

void string_copy(struct string_t *dest, struct string_t *src);

void string_write(struct string_t *string, char *text, unsigned int msgSize);

bool string_cut(struct string_t *string, unsigned int cutFromStart, unsigned int cutFromEnd);

void string_destroy(struct string_t *string);

void string_print(struct string_t *string);

struct buffer_queue_t *buffer_queue_create(enum buffer_queue_mode_t mode);

void buffer_queue_destroy(struct buffer_queue_t *bufferQueue);

bool buffer_queue_write(struct buffer_queue_t *bufferQueue, char *message, unsigned int msgSize);

bool buffer_queue_read(struct buffer_queue_t *bufferQueue, struct string_t *textRead);

unsigned int buffer_queue_get_read_size(struct buffer_queue_t *bufferQueue);

bool buffer_queue_reset(struct buffer_queue_t *bufferQueue);

struct buffer_queue_t *buffer_queue_split(struct string_t *string, char pattern);

void buffer_queue_concatenate(struct buffer_queue_t *bufferQueue, struct string_t *string);

void buffer_queue_print(struct buffer_queue_t *bufferQueue);


#endif