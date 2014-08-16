#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct cluster_s cluster_t;
typedef struct cluster_node_s cluster_node_t;
typedef struct proximity_node_s proximity_node_t;

struct cluster_s {
     int num_items;
     int num_clusters;
     cluster_node_t *nodes;
};

struct cluster_node_s {
     int is_used;
     int is_root;
     int height;
     char *label;
     int num_items;
     cluster_node_t **items;
     cluster_node_t **merged;
     proximity_node_t *proximity;
};

struct proximity_node_s {
     int target;
     float distance;
     proximity_node_t *next, *prev;
};

cluster_node_t *add_leaf(cluster_t *cluster, const char *label) {
     cluster_node_t *t = &cluster->nodes[cluster->num_items];
     if ((t->label = (char *)
          calloc(strlen(label), sizeof(char))) == NULL) {
          fprintf(stderr, "Failed to allocate label.\n");
          t = NULL;
     } else {
          strlcpy(t->label, label, sizeof(label));
          t->is_used = 1;
          t->is_root = 1;
          cluster->num_items++;
          cluster->num_clusters++;
     }
     return t;
}

void free_proximity_list(proximity_node_t *node) {
     proximity_node_t *t;
     while (node) {
          t = node->next;
          free(node);
          node = t;
     }
}

void free_cluster(cluster_t * cluster) {
     int i;
     cluster_node_t node;
     for (i = 0; i < cluster->num_items; ++i) {
          node = cluster->nodes[i];
          if (node.label != NULL)
               free(node.label);
          if (node.proximity != NULL)
               free_proximity_list(node.proximity);
     }
     free(cluster);
}

void insert_before(proximity_node_t *current, proximity_node_t *proximity,
                   cluster_node_t *node) {
     proximity->next = current;
     if (current->prev == NULL) {
          node->proximity = proximity;
     } else {
          current->prev->next = proximity;
          proximity->prev = current->prev;
     }
     current->prev = proximity;
}

void insert_after(proximity_node_t *current, proximity_node_t *proximity) {
     proximity->prev = current;
     current->next = proximity;
}

void insert_sorted(cluster_node_t *node, proximity_node_t *proximity) {
     proximity_node_t *temp = node->proximity;
     while(temp->next != NULL) {
          if (temp->distance >= proximity->distance) {
               insert_before(temp, proximity, node);
               return;
          }
          temp = temp->next;
     }
     if (proximity->distance < temp->distance)
          insert_before(temp, proximity, node);
     else
          insert_after(temp, proximity);
}

proximity_node_t *add_proximity(cluster_node_t *node, int index, int target,
                                const float **distance) {
     proximity_node_t *proximity;
     if ((proximity = (proximity_node_t *)
          calloc(1, sizeof(proximity_node_t))) == NULL) {
          fprintf(stderr, "Failed to allocate proximity node.\n");
     } else {
          proximity->target = target;
          proximity->distance = distance[index][target];
          if (node->proximity == NULL)
               node->proximity = proximity;
          else
               insert_sorted(node, proximity);
     }
     return proximity;
}

cluster_t *update_proximity(cluster_t *cluster, int index,
                            const float **distance) {
     int i = 1, target = index;
     cluster_node_t node = cluster->nodes[index], temp;
     if (node.is_used) {
          while (i++ < cluster->num_clusters) {
               temp = cluster->nodes[--target];
               if (!temp.is_used) {
                    fprintf(stderr, "Empty node found at %d.\n", index);
                    cluster = NULL;
                    break;
               }
               if (temp.is_root)
                    add_proximity(&node, index, target, distance);
          }
     } else {
          fprintf(stderr, "Supplied node with index %d is invalid\n", index);
          cluster = NULL;
     }
     return cluster;
}

cluster_t *add_leaves(cluster_t *cluster, const char *labels[],
                      const float **distance) {
     int i = 0;
     cluster_node_t *node;
     for (i = 0; i < cluster->num_items; ++i) {
          node = add_leaf(cluster, labels[i]);
          if (node == NULL) {
               free_cluster(cluster);
               cluster = NULL;
               break;
          } else
               update_proximity(cluster, i, distance);
     }
     return cluster;
}

cluster_t *merge_clusters(cluster_t *cluster) {
     return cluster;
}

cluster_t *agglomerate(int num_items, const char *labels[], const float **distance) {
     cluster_t *cluster;
     if ((cluster = (cluster_t *) calloc(1, sizeof(cluster_t))) == NULL)
          fprintf(stderr, "Failed to allocate cluster.\n");
     else {
          cluster->num_items = num_items;
          cluster->num_clusters = 0;
          if ((cluster->nodes = (cluster_node_t *)
               calloc(2 * num_items - 1, sizeof(cluster_node_t))) == NULL) {
               fprintf(stderr, "Failed to allocate cluster nodes.\n");
               free(cluster);
               cluster = NULL;
          } else {
               add_leaves(cluster, labels, distance);
               merge_clusters(cluster);
          }
     }
     return cluster;
}

int main(int argc, char **argv) {
     return 0;
}
