#include <stdio.h>
#include <stdlib.h>
#include "list.h"

// Check if the list L is empty
int is_empty(List L) { return !L; }
 
// Print the List L
void print(List L) { 
  	List p;  	p=L;
  	printf("Indexes : ");  	while (!is_empty(p)) {
  		printf("%d ;", p->index);
  		p=p->next;
  }
  printf("\n");
}

// Add a element with the value number 
List append(int number, List End){
	List p=calloc(1,sizeof(*p));
	if (p==NULL) return NULL;
	
	if (End == NULL){
		End = calloc(1,sizeof(*End));
		End->index = number;
		return End;
	}
	
	p->index = number;
	p->next=NULL;
	End->next=p;
	return End;
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
	