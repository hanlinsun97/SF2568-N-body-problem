#ifndef _LIST
#define _LIST

struct link {	int number;    struct link *next;} ;
typedef struct link* List;

int is_empty(List L);
void print(List L);
List append(int number, List L);

#endif
