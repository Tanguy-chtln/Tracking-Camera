#include <pthread.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include "SerialQueue.h"


struct string_chain_t {
    struct string_chain_t *next;
    struct string_chain_t *last;
    struct string_t string;
};

struct buffer_queue_t {
    struct string_chain_t *readPtr;
    struct string_chain_t *first;
    struct string_chain_t *last;
    enum buffer_queue_mode_t mode;
    pthread_mutex_t mutex;
    unsigned int readSize;
};


void string_init(struct string_t *string) {
    string->text = NULL;
    string->textLength = 0;
}

void string_copy(struct string_t *dest, struct string_t *src) {
    if (dest->text != NULL)
        free(dest->text);
    dest->text = calloc((src->textLength+1), sizeof(char));
    dest->textLength = src->textLength;
    memcpy(dest->text, src->text, src->textLength);
}

void string_write(struct string_t *string, char *text, unsigned int msgSize) {
    if (msgSize == 0)
        return;
    if (string->text != NULL)
        free(string->text);
    
    string->text = malloc(msgSize+1);
    memcpy(string->text, text, msgSize);
    string->textLength = msgSize;
    string->text[msgSize] = '\0';
}

bool string_cut(struct string_t *string, unsigned int cutFromStart, unsigned int cutFromEnd) {
    if (string->textLength == cutFromStart+cutFromEnd){
        free(string->text);
        string->text = calloc(1, sizeof(char));
        string->textLength = 0;
    }

    if (string->textLength <= cutFromStart+cutFromEnd)
        return false;
    char *textTmp = calloc(string->textLength-(cutFromStart+cutFromEnd), sizeof(char));
    memcpy(textTmp, string->text+cutFromStart, string->textLength-(cutFromStart+cutFromEnd));
    string_write(string, textTmp, string->textLength-(cutFromStart+cutFromEnd));
    free(textTmp);
    return true;
}

void string_destroy(struct string_t *string) {
    if (string->text != NULL) {
        free(string->text);
        string->text = NULL;
        string->textLength = 0;
    }
}

bool sublist_idx(struct string_t *string, unsigned int stringOffset, char **subListStartIdx, char **subListEndIdx) {
    if (stringOffset >= string->textLength)
        return false;
    char startPattern = '[';
    char endPattern = ']';
    *subListStartIdx = memchr(string->text + stringOffset, startPattern, string->textLength-stringOffset);
    if (*subListStartIdx == NULL)
        return false;

    *subListEndIdx = memchr(*subListStartIdx, endPattern, string->textLength - (unsigned int) (*subListStartIdx - string->text));

    if (
        *subListEndIdx != NULL && 
        ((*subListStartIdx > string->text && *(*subListStartIdx-1) == ' ') || *subListStartIdx == string->text) &&
        ((*subListEndIdx < string->text + string->textLength-1  && ((*subListEndIdx)[1] == ' ' || (*subListEndIdx)[1] == '\n' || (*subListEndIdx)[1] == '\r')) || *subListEndIdx == string->text+string->textLength-1 )
        ) 
    {
        return true;
    }
    return false;

}


unsigned int string_count(struct string_t *string, char pattern, char **idxBuff) {
    unsigned int nAnsw = 0;
    unsigned int lengthRead = 0;
    char *searchPos;
    char *subListStartIdx;
    char *subListEndIdx;

    searchPos = (char *) memchr(string->text + lengthRead, pattern, string->textLength - lengthRead);
    if (sublist_idx(string, lengthRead, &subListStartIdx, &subListEndIdx) && subListStartIdx <= searchPos) {
        lengthRead = (subListEndIdx - string->text) + 1;
        if (subListEndIdx < string->text+string->textLength-1) {
            searchPos = subListEndIdx + 1;
            lengthRead = (subListEndIdx - string->text) + 2;
        } else {
            searchPos = subListEndIdx ;
            lengthRead = (subListEndIdx - string->text)+1;
        }
    }
    else
        lengthRead = (searchPos - string->text) + 1;
    while (lengthRead < string->textLength && searchPos != NULL) {
        idxBuff[nAnsw++] = searchPos;
        searchPos = (char *) memchr(string->text + lengthRead, pattern, string->textLength - lengthRead);
        if (sublist_idx(string, lengthRead, &subListStartIdx, &subListEndIdx) && subListStartIdx < searchPos) {
            if (subListEndIdx < string->text+string->textLength-1){
                searchPos = subListEndIdx + 1;
                lengthRead = (subListEndIdx - string->text) + 2;
            } else {
                searchPos = subListEndIdx;
                lengthRead = (subListEndIdx - string->text)+1;
            }
        }
        else
            lengthRead = (searchPos - string->text) + 1;
    }
    if (lengthRead < string->textLength && searchPos != NULL) 
        idxBuff[nAnsw++] = searchPos;
    
    return nAnsw;
}

void string_print(struct string_t *string) {
    char *printBuff;
    printBuff = malloc(string->textLength+1);
    memcpy(printBuff, string->text, string->textLength);
    printBuff[string->textLength] = 0;
    printf("%s\n", printBuff);
    free(printBuff);
}



struct buffer_queue_t *buffer_queue_create(enum buffer_queue_mode_t mode) {
    if (mode != QUEUE && mode != KEEP_DATA)
        return NULL;
    struct buffer_queue_t *bufferQueue = calloc(1, sizeof(struct buffer_queue_t));
    bufferQueue->readPtr = NULL;
    bufferQueue->first = NULL;
    bufferQueue->last = NULL;
    bufferQueue->mode = mode;
    bufferQueue->readSize = 0;
    pthread_mutex_init(&bufferQueue->mutex, NULL);
    return bufferQueue;
}

void buffer_queue_destroy(struct buffer_queue_t *bufferQueue) {
    struct string_chain_t *chain = bufferQueue->first;
    while ( chain != bufferQueue->last) {
        chain = bufferQueue->first->next;
        string_destroy(&bufferQueue->first->string);
        free(bufferQueue->first);
        bufferQueue->first = chain;
    }
    if (bufferQueue->first != NULL) {
        string_destroy(&bufferQueue->first->string);
        free(bufferQueue->first);
    }
    pthread_mutex_destroy(&bufferQueue->mutex);
    free(bufferQueue);
}

bool buffer_queue_write(struct buffer_queue_t *bufferQueue, char *message, unsigned int msgSize) {
    bool answ;
    pthread_mutex_lock(&bufferQueue->mutex);
    if (msgSize == 0)
        answ = false;
    answ = true;
    struct string_chain_t *nextNode = calloc(1, sizeof(struct string_chain_t));
    string_write(&nextNode->string, message, msgSize);
    
    if (bufferQueue->last == NULL) {
        bufferQueue->last = nextNode;
        bufferQueue->first = nextNode;
        bufferQueue->readPtr = nextNode;
    } else {
        bufferQueue->last->next = nextNode;
        bufferQueue->last = nextNode;
    }

    bufferQueue->readSize++;
    pthread_mutex_unlock(&bufferQueue->mutex);
    return answ;
}


bool buffer_queue_read(struct buffer_queue_t *bufferQueue, struct string_t *textRead) {
    bool answ;
    pthread_mutex_lock(&bufferQueue->mutex);
    if (bufferQueue->readPtr != NULL) {
        answ = true;
        string_copy(textRead, &bufferQueue->readPtr->string);
        
        bufferQueue->readPtr = bufferQueue->readPtr->next;
        if (bufferQueue->mode != KEEP_DATA) {
            string_destroy(&bufferQueue->first->string);
            free(bufferQueue->first);
            if (bufferQueue->first != bufferQueue->last)
                bufferQueue->first = bufferQueue->readPtr;
            else {
                bufferQueue->first = bufferQueue->last = bufferQueue->readPtr = NULL;
            }

        }

        bufferQueue->readSize--;
    }
    else {
        answ = false;
    }
    pthread_mutex_unlock(&bufferQueue->mutex);
    return answ;
}


unsigned int buffer_queue_get_read_size(struct buffer_queue_t *bufferQueue) {
    return bufferQueue->readSize;
}


bool buffer_queue_reset(struct buffer_queue_t *bufferQueue) {
    if (bufferQueue->mode != KEEP_DATA) 
        return false;
    
    pthread_mutex_lock(&bufferQueue->mutex);
    bufferQueue->readPtr = bufferQueue->first;
    pthread_mutex_unlock(&bufferQueue->mutex);
    return true;
}





struct buffer_queue_t *buffer_queue_split(struct string_t *string, char pattern) {
    struct buffer_queue_t *bufferQueueSplit = buffer_queue_create(KEEP_DATA);
    if (string->textLength == 0)
        return bufferQueueSplit;

    char **splitPoints = malloc(sizeof(char *) * string->textLength);
    unsigned int nbSplits = string_count(string, pattern, splitPoints);
    

    char *textRead = string->text;
    
    
    for (unsigned int i=0; i<nbSplits; i++) {
        buffer_queue_write(bufferQueueSplit, textRead, (unsigned int) (splitPoints[i] - textRead));
        textRead = splitPoints[i] + 1;
    }
    buffer_queue_write(bufferQueueSplit, textRead, (unsigned int) ((string->text + string->textLength) - textRead));
    free(splitPoints);
    return bufferQueueSplit;
}

void buffer_queue_concatenate(struct buffer_queue_t *bufferQueue, struct string_t *string) {
    pthread_mutex_lock(&bufferQueue->mutex);
    struct string_chain_t *ptr = bufferQueue->readPtr;
    char *stringPtr;
    free(string->text);
    string->textLength = 0;
    if (ptr == NULL) {
        string->text = calloc(1, sizeof(char));
        return;
    }

    while(ptr != bufferQueue->last) {
        string->textLength += ptr->string.textLength;
        ptr = ptr->next;
    }
    string->textLength += ptr->string.textLength;
    ptr = bufferQueue->readPtr;
    string->text = malloc(string->textLength+1);
    stringPtr = string->text;

    while(ptr != bufferQueue->last) {
        memcpy(stringPtr, ptr->string.text, ptr->string.textLength);
        stringPtr += ptr->string.textLength;
        ptr = ptr->next;
    }
    memcpy(stringPtr, ptr->string.text, ptr->string.textLength);
    stringPtr += ptr->string.textLength;

    *stringPtr = '0';
    pthread_mutex_unlock(&bufferQueue->mutex);
}


void buffer_queue_print(struct buffer_queue_t *bufferQueue) {
    struct string_chain_t * printPtr=bufferQueue->readPtr;
    while (printPtr != NULL) {
        string_print(&printPtr->string);
        printPtr = printPtr->next;
    }
}

