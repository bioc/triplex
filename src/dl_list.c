/**
 * Triplex package
 * Result export and filtering
 *
 * @author  Tomas Martinek
 * @date    2012/10/15
 * @file    dl_list.c
 * @package triplex
 */

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dl_list.h"

/* #define DEBUG */

/*************************************************************************************************/
/*************************************************************************************************/
int test_duplication(t_dl_node *lst, t_dl_data *new) {
   
   /* the same positions and scores */
   if ((lst->data.start == new->start) && (lst->data.lstart == new->lstart) &&
      (lst->data.lend == new->lend) && (lst->data.end == new->end) &&
      (lst->data.score == new->score)) {
#ifdef DEBUG
      Rprintf("Duplication (the same score)\n");
#endif
      return 1;
   }
      
   /* score differs - if new one is better then update  */
   if ((lst->data.start == new->start) && (lst->data.end == new->end)) {
      if (lst->data.score < new->score)
         lst->data = *new;
#ifdef DEBUG
      Rprintf("Duplication (lower score)\n");
#endif
      return 1;
   }


   return 0;
}

/*************************************************************************************************/
/*************************************************************************************************/
int test_included(t_dl_list *list, t_dl_node *node, t_dl_data *new) {

   t_dl_node *temp = node, *pointer;
   
   /* next node is the first one with possible larger range */
   pointer = node->next;

   /* iterate through the nodes with the same start position and larger end 
   (to right) */
   while((pointer != NULL) && (pointer->data.start == new->start)) {
      if(pointer->data.score >= new->score) {
#ifdef DEBUG
         Rprintf("Node is included (the same start)\n");
#endif
        return 1;
      }
      pointer = pointer->next;
   }

   pointer = temp;

   /* find the first node with lower start position (to left) */
   while(pointer->data.start == new->start)
      pointer = pointer->prev;

   /* iterate through the nodes with lower start position (to the left) */
   while(pointer->data.start >= new->start-list->max_len) {
      if((pointer->data.end >= new->end) &&
         (pointer->data.score >= new->score)) {
#ifdef DEBUG
         Rprintf("Node is included (lower start)\n");
#endif
         return 1;
      }
      pointer = pointer->prev;
   }

   return 0;
}

/*************************************************************************************************/
/*************************************************************************************************/
int test_include(t_dl_list *list, t_dl_node *node) {

   t_dl_data *new = &node->data;
   t_dl_node *temp, *pointer = node->prev;

   /* iterate through the nodes with the same start position and lower end
   (to left) */
   while(pointer->data.start == new->start) {
      temp = pointer->prev;
      if(pointer->data.score <= new->score) {
#ifdef DEBUG
         Rprintf("Node includes node (the same start)\n");
#endif
         dl_list_delete(list, pointer);
      }
      pointer = temp;
   }

   pointer = node->next;

   /* find the first node with higher start position (to right) */
   while((pointer != NULL) && (pointer->data.start == new->start))
      pointer = pointer->next;

   /* iterate through the nodes with higher start position and lower end (to
   right) */
   while((pointer != NULL) && (pointer->data.end <= new->end)){
      temp = pointer->next;
      if (pointer->data.score <= new->score) {
#ifdef DEBUG
         Rprintf("Node includes node (higher start)\n");
#endif
         dl_list_delete(list, pointer);
      }
      pointer = temp;
   }

   return 0;
}

/*************************************************************************************************/
/*************************************************************************************************/
int dl_list_insert(t_dl_list *list, t_dl_data data)
{
        t_dl_node *temp, *pointer = list->last;

        /* Iterate through the list till we encounter right position */
        while((data.start < pointer->data.start) || 
             ((data.start == pointer->data.start) && (data.end < pointer->data.end))) {
                pointer = pointer->prev;
        }
        /* Do not insert the same data */
        if(test_duplication (pointer, &data))
           return 0;
        if(test_included(list, pointer, &data))
           return 0;

        /* Allocate memory for the new node and put data in it.*/
        temp = pointer->next;
        pointer->next = (t_dl_node *)malloc(sizeof(t_dl_node));
        (pointer->next)->prev = pointer;
        pointer = pointer->next;
        pointer->data = data;
        pointer->next = temp;
        if(temp != NULL)
           temp->prev=pointer;
        if(pointer->prev == list->last)
           list->last = pointer;
        list->size++;

        /* Remove all included nodes with lower score */
        test_include(list, pointer);

        return 0;
}

/*************************************************************************************************/
/*************************************************************************************************/
int dl_node_overlap(t_dl_node *n1, t_dl_node *n2)
{
   int overlap, whole;


   if (n1->data.end > n2->data.start) {
      overlap = n1->data.end - n2->data.start; 
      whole = n2->data.end - n1->data.start;

      if (((float)overlap/(float)whole) >= 0.8)
         return 1;
   }

   return 0;
}

/*************************************************************************************************/
/*************************************************************************************************/
void local_group_filter(t_dl_list *list, t_dl_node *start, t_dl_node *end)
{
   t_dl_node *pointer, *temp;
   t_dl_node *new_start, *new_end;
   int change;

#ifdef DEBUG
   Rprintf("Group: (%d,%d) - (%d,%d)\n", start->data.start, start->data.end,
                                        end->data.start, end->data.end);
#endif

   do {
      /* iterate through neighbour elements and mark those for deletion */
      pointer = start;
      change = 0;
      while (pointer != end) {
         if (dl_node_overlap(pointer, pointer->next) 
            /* && (pointer->data.score != (pointer->next)->data.score) */
            ) {
            /* mark first element */
            if(pointer->data.score < (pointer->next)->data.score)
               pointer->data.type = -1;
            /* mark second element */
            else
               (pointer->next)->data.type = -1;
            change = 1;
         }
         pointer = pointer->next;
      }

      if (change) {
         /* iterate through group elements and delete marked ones */
         pointer = start;
         new_start = start;
         new_end = end;
         while (pointer != end->next) {
#ifdef DEBUG
            Rprintf("Element: (%d,%d) - %d\n", pointer->data.start,
               pointer->data.end, pointer->data.type);
#endif
            temp = pointer;
            pointer = pointer->next;
            if (temp->data.type == -1) {
               if(new_start == temp)
                  new_start = temp->next;
               if(new_end == temp)
                  new_end = temp->prev;
               dl_list_delete(list, temp);
            }
         }
         start = new_start;
         end = new_end;
      }

   }while(change);
}

/*************************************************************************************************/
/*************************************************************************************************/
void dl_list_group_filter(t_dl_list *list) 
{
   t_dl_node *pointer = (list->first)->next;
   t_dl_node *group_start, *group_end;

   while(pointer!=NULL) {

      /* Group detection */
      group_start = pointer;
      while((pointer->next!=NULL) && dl_node_overlap(pointer, pointer->next))
         pointer = pointer->next;
      group_end = pointer;

      pointer = pointer->next;
         
      /* Group filtration */
      if (group_start != group_end)
         local_group_filter(list, group_start, group_end);
   }
}

/*************************************************************************************************/
/*************************************************************************************************/
void dl_list_get_head_range(t_dl_list *list, int *start, int *end) 
{
   t_dl_node *pointer = (list->first)->next;

   if (pointer != NULL) {
      *start = pointer->data.start;
      *end = pointer->data.end;
   }
}

/*************************************************************************************************/
/*************************************************************************************************/
void dl_list_merge_sort(t_dl_list *list_arr, t_dl_list *list_out, int num)
{
   int i, min_list, min_start, min_end, act_start, act_end, items = 0;
   t_dl_node *pointer;

   /* Initialization of output list */
   dl_list_init(list_out, list_arr[0].max_len);

   /* Overall number of items */
   for(i=0; i<num; i++)
      items += list_arr[i].size;

#ifdef DEBUG
   Rprintf("Items %d\n", items);
#endif

   while(items > 0) {
      /* Find a list with the lowest element */
      min_list = -1;
      for(i=0; i<num; i++) {
         if(list_arr[i].size!=0) {
            if (min_list == -1) {
               dl_list_get_head_range(&list_arr[i], &min_start, &min_end);
               min_list = i;
            } else {
               dl_list_get_head_range(&list_arr[i], &act_start, &act_end);
               if ((act_start < min_start) || 
                   ((act_start == min_start) && (act_end < min_end))){
                  min_start = act_start;
                  min_end = act_end;
                  min_list = i;
               }
            }
         }
      }

#ifdef DEBUG
      Rprintf("Move element form list %d\n", min_list);
#endif

      /* Move the lowest element from input list to the end of output list */
      pointer = (list_arr[min_list].first)->next;
      (pointer->prev)->next = pointer->next;
      list_arr[min_list].size--;
      if (pointer->next != NULL)
         (pointer->next)->prev = pointer->prev;
      else
         list_arr[min_list].last = pointer->prev;

      (list_out->last)->next = pointer;
      pointer->prev = list_out->last;
      pointer->next = NULL;
      list_out->last = pointer;
      list_out->size++;

      items--;
   }
}

/*************************************************************************************************/
/*************************************************************************************************/
void dl_list_delete(t_dl_list *list, t_dl_node *node)
{
        t_dl_node *pointer = node->prev;

        pointer->next = node->next;

        if (node->next != NULL)
           (node->next)->prev = pointer;
        else
           list->last = pointer;

        free(node);
        list->size--;
}

/*************************************************************************************************/
/*************************************************************************************************/
void dl_list_init(t_dl_list *list, int max_len)
{
        list->first = (t_dl_node *)malloc(sizeof(t_dl_node));
        (list->first)->data.start = -2*max_len-1;
        (list->first)->data.end = -1;
        (list->first)->data.score = 0;
        (list->first)->next = NULL;
        (list->first)->prev = NULL;
        list->last = list->first;
        list->size = 0;
        list->max_len = max_len;
}

/*************************************************************************************************/
/*************************************************************************************************/
void dl_list_free(t_dl_list *list)
{
   t_dl_node *temp, *pointer = list->first;

   while(pointer!=NULL) {
       temp = pointer;
      pointer = pointer->next;
       free(temp);
   }
   list->first = NULL;
   list->last = NULL;
   list->size = 0;
   list->max_len = 0;
}

/*************************************************************************************************/
/*************************************************************************************************/
void dl_list_print(t_dl_list *list)
{
        t_dl_node *pointer = (list->first)->next;
        /* t_dl_node *pointer = list->first; */

        while(pointer!=NULL) {
            Rprintf("[%d,%d:%d]", pointer->data.start, 
                                 pointer->data.end, 
                                 pointer->data.score);
/*
            Rprintf("(");
            if (pointer->prev != NULL)
               Rprintf("[%d,%d:%d]", (pointer->prev)->data.start, 
                                    (pointer->prev)->data.end,
                                    (pointer->prev)->data.score);
            else
               Rprintf("NULL");
            Rprintf(",");
            if (pointer->next != NULL)
               Rprintf("[%d,%d:%d]", (pointer->next)->data.start, 
                                    (pointer->next)->data.end,
                                    (pointer->next)->data.score);
            else
               Rprintf("NULL");
            Rprintf(") ");
*/
            pointer = pointer->next;
        }
        Rprintf("\nList size: %d", list->size);
        Rprintf("\n");
}

/*************************************************************************************************/
/*************************************************************************************************/
void assign_data(t_dl_data *data, int start, int lstart, 
                 int lend, int end, int score) {
   data->start = start;
   data->lstart = lstart;
   data->lend = lend;
   data->end = end;
   data->score = score;
   data->type = 0;
}

/*************************************************************************************************/
/*************************************************************************************************/
/*int main() {
        t_dl_data data;
        t_dl_list list;
        t_dl_list list_arr[4];
        int i,j;

        dl_list_init(&list, 15);*/

        /* Duplication test */
/*
        assign_data(&data, 10, 15, 20, 25, 20);
        insert(&list,data);
        print(&list);
        assign_data(&data, 10, 15, 20, 25, 20);
        insert(&list,data);
        print(&list);
        assign_data(&data, 10, 16, 19, 25, 25);
        insert(&list,data);
        print(&list);
*/
        /* Included test */
/*
        assign_data(&data, 10, 15, 20, 25, 20);
        insert(&list,data);
        print(&list);
        assign_data(&data, 10, 15, 20, 24, 15);
        insert(&list,data);
        print(&list);
        assign_data(&data, 11, 15, 20, 24, 15);
        insert(&list,data);
        print(&list);
        assign_data(&data, 11, 15, 20, 24, 25);
        insert(&list,data);
        print(&list);
        assign_data(&data, 10, 15, 20, 23, 25);
        insert(&list,data);
        print(&list);
*/
        /* Include test */
/*
        assign_data(&data, 10, 15, 20, 25, 20);
        insert(&list,data);
        print(&list);
        assign_data(&data, 12, 15, 20, 27, 20);
        insert(&list,data);
        print(&list);
        assign_data(&data, 15, 17, 20, 30, 20);
        insert(&list,data);
        print(&list);
        assign_data(&data, 10, 15, 20, 30, 20);
        insert(&list,data);
        print(&list);
*/

        /* Random test */
/*
        srand(time(NULL));
        for(i=0; i<20000; i++) {
            data.start = rand() % 1000;
            data.lstart = data.start + rand() % 20;
            data.lend = data.lstart + rand() % 20;
            data.end = data.lend + rand() % 20;
            data.score = 10 + rand() % 20;
            dl_list_insert(&list,data);
        }
        dl_list_print(&list);
        dl_list_free(&list);
*/

        /* List merge test */
/*
        srand(time(NULL));
        Rprintf("============================================\n");
        for (i=0; i<4; i++) {
           dl_list_init(&list_arr[i], 15);

           for(j=0; j<4; j++) {
               data.start = rand() % 1000;
               data.lstart = data.start + rand() % 20;
               data.lend = data.lstart + rand() % 20;
               data.end = data.lend + rand() % 20;
               data.score = 10 + rand() % 20;
               dl_list_insert(&list_arr[i],data);
           }

           dl_list_print(&list_arr[i]);
        }
        Rprintf("============================================\n");

        dl_list_merge_sort(list_arr, &list, 4);
        dl_list_print(&list);
        dl_list_free(&list);
*/

        /* Group filtration test - Simple */
/*
        assign_data(&data, 10, 15, 20, 45, 22);
        dl_list_insert(&list, data);
        assign_data(&data, 11, 15, 20, 46, 21);
        dl_list_insert(&list, data);
        assign_data(&data, 12, 15, 20, 47, 20);
        dl_list_insert(&list, data);
        assign_data(&data, 13, 15, 20, 48, 21);
        dl_list_insert(&list, data);
        dl_list_print(&list);
        dl_list_group_filter(&list);
        dl_list_print(&list);
        dl_list_free(&list);
*/
        /* Group filtration test - Random */
/*        srand(time(NULL));
        for(i=0; i<20000; i++) {
            data.start = rand() % 1000;
            data.lstart = data.start + rand() % 50;
            data.lend = data.lstart + rand() % 50;
            data.end = data.lend + rand() % 50;
            data.score = 10 + rand() % 50;
            dl_list_insert(&list,data);
        }
        dl_list_print(&list);
        dl_list_group_filter(&list);
        dl_list_print(&list);
        dl_list_free(&list);

        return 0;
}*/
