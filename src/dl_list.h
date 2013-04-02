/**
 * Triplex package
 * Header file for result export and filtering
 *
 * @author  Tomas Martinek
 * @date    2012/10/15
 * @file    dl_list.h
 * @package triplex
 */

#ifndef DL_LIST_H
#define DL_LIST_H

typedef struct DL_Data {
	int type;
	int start;
	int end;
	int lstart;
	int lend;
	int score;
	double pvalue;
	int insdel;
	int strand;
} t_dl_data;

typedef struct DL_Node {
        struct DL_Data data;
        struct DL_Node *next;
        struct DL_Node *prev;
}t_dl_node;

typedef struct {
        int    size;
        int    max_len;
        struct DL_Node *first;
        struct DL_Node *last;
}t_dl_list;

void dl_list_init(t_dl_list *list, int max_len);
void dl_list_delete(t_dl_list *list, t_dl_node *node);
int dl_list_insert(t_dl_list *list, t_dl_data data);
void dl_list_free(t_dl_list *list);
void dl_list_merge_sort(t_dl_list *list_arr, t_dl_list *list_out, int num);
void dl_list_group_filter(t_dl_list *list);

#endif // DL_LIST_H
