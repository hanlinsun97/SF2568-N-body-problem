#include <stdio.h>
#include <stdlib.h>
#include "list.h"

int is_empty(List L) { return !L; }
 
void print(List L) { 
  	List p;  	p=L;  	while (!is_empty(p)) {
  		printf("%d \n", p->number);
  		p=p->next;
  }
}

List append(int number, List L){
	List c;
	List p=calloc(1,sizeof(*p));
	if (p==NULL) return NULL;
	
	for(c=L; !is_empty(c->next);c=c->next);
	
	p->number = number;
	p->next=NULL;
	c->next=p;
	return L;
}

// TO TEST THE LISTS
//int main(){
//	List test = calloc(1,sizeof(*test));
//	test->next = NULL;
//	test = append(1,test);
//	test = append(4,test);
//	test = append(3,test);
//
//	print(test);
//	return 0;
//}
	