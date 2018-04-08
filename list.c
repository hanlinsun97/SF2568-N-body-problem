#include <stdio.h>
#include <stdlib.h>
#include "list.h"

// Check if the list L is empty
int is_empty(List L) { return !L; }
 
// Print the List L
void print(List L) { 
  	List p;  	p=L;  	while (!is_empty(p)) {
  		printf("%d \n", p->index);
  		p=p->next;
  }
}

// Add a element with the value number at the end of the List L
List append(int number, List L){
	List c;
	List p=calloc(1,sizeof(*p));
	if (p==NULL) return NULL;
	
	if (L == NULL){
		L = calloc(1,sizeof(*L));
		L->index = number;
		return L;
	}

	for(c=L; !is_empty(c->next);c=c->next);
	
	p->index = number;
	p->next=NULL;
	c->next=p;
	return L;
}

// TO TEST THE LISTS
//int main(){
//	List test = NULL;
//	test = append(1,test);
//	test = append(4,test);
//	test = append(3,test);
//
//	print(test);
//	return 0;
//}
	