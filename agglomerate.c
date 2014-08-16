#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct cluster_s cluster_t;
typedef struct cluster_node_s cluster_node_t;
typedef struct proximity_node_s proximity_node_t;
typedef struct item_s item_t;

struct cluster_s {
     int num_items;
     int num_clusters;
     cluster_node_t *nodes;
     float **distances;
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

struct item_s {
     float x, y;
     char *label;
};

cluster_node_t *add_leaf(cluster_t *cluster, const char *label) {
     int retval = 0;
     cluster_node_t *t = &(cluster->nodes[cluster->num_clusters]);
     t->label = (char *) calloc(strlen(label), sizeof(char));
     if (t->label) {
          strlcpy(t->label, label, sizeof(label));
          t->is_used = 1;
          t->is_root = 1;
          cluster->num_clusters++;
     } else {
          fprintf(stderr, "Failed to allocate label.\n");
          t = NULL;
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
     if (cluster) {
          if (cluster->nodes) {
               for (i = 0; i < cluster->num_items; ++i) {
                    node = cluster->nodes[i];
                    if (node.label)
                         free(node.label);
                    if (node.proximity)
                         free_proximity_list(node.proximity);
               }
               free(cluster->nodes);
          }
          if (cluster->distances)
               free(cluster->distances);
          free(cluster);
     }
}

void insert_before(proximity_node_t *current, proximity_node_t *proximity,
                   cluster_node_t *node) {
     proximity->next = current;
     if (current->prev) {
          current->prev->next = proximity;
          proximity->prev = current->prev;
     } else
          node->proximity = proximity;
     current->prev = proximity;
}

void insert_after(proximity_node_t *current, proximity_node_t *proximity) {
     proximity->prev = current;
     current->next = proximity;
}

void insert_sorted(cluster_node_t *node, proximity_node_t *proximity) {
     proximity_node_t *temp = node->proximity;
     while (temp->next) {
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

proximity_node_t *add_proximity(cluster_t *cluster, int index, int target) {
     proximity_node_t *proximity = (proximity_node_t *)
          calloc(1, sizeof(proximity_node_t));
     if (proximity) {
          proximity->target = target;
          proximity->distance = cluster->distances[index][target];
          cluster_node_t *node = &(cluster->nodes[index]);
          if (node->proximity)
               insert_sorted(node, proximity);
          else
               node->proximity = proximity;
     } else 
          fprintf(stderr, "Failed to allocate proximity node.\n");
     return proximity;
}

cluster_t *update_proximity(cluster_t *cluster, int index) {
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
                    add_proximity(cluster, index, target);
          }
     } else {
          fprintf(stderr, "Supplied node with index %d is invalid\n", index);
          cluster = NULL;
     }
     return cluster;
}

float euclidean_distance(item_t a, item_t b) {
     return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
}

void fill_eucliden_distances(float **matrix, int num_items,
                             const item_t items[]) {
     int i, j;
     for (i = 0; i < num_items; ++i)
          for (j = 0; j < num_items; ++j) {
               matrix[i][j] = euclidean_distance(items[i], items[j]);
               matrix[j][i] = matrix[i][j];
          }
}

float **generate_distance_matrix(int num_items, const item_t items[]) {
     int i, j;
     float **matrix = matrix = (float **) calloc(num_items, sizeof(float *));
     if (matrix) {
          for (i = 0; i < num_items; ++i) {
               matrix[i] = (float *) calloc(num_items, sizeof(float));
               if (!matrix[i]) {
                    fprintf(stderr, "Failed to allocate distance matrix.\n");
                    num_items = i;
                    for (i = 0; i < num_items; ++i)
                         free(matrix[i]);
                    free(matrix);
                    matrix = NULL;
                    break;
               }    
          }
          if (matrix)
               fill_eucliden_distances(matrix, num_items, items);
     } else
          fprintf(stderr, "Failed to allocate distance matrix.\n");
     return matrix;
}

cluster_t *add_leaves(cluster_t *cluster, const item_t items[]) {
     int i = 0;
     cluster_node_t *node;
     for (i = 0; i < cluster->num_items; ++i) {
          if (add_leaf(cluster, items[i].label))
               update_proximity(cluster, i);
          else {
               cluster = NULL;
               break;
          }
     }
     return cluster;
}

cluster_t *merge_clusters(cluster_t *cluster) {
     
     return cluster;
}

cluster_t *agglomerate(int num_items, const item_t items[]) {
     cluster_t *cluster = (cluster_t *) calloc(1, sizeof(cluster_t));
     if (cluster) {
          cluster->nodes = (cluster_node_t *)
               calloc(2 * num_items - 1, sizeof(cluster_node_t));
          if (cluster->nodes) {
               cluster->distances = generate_distance_matrix(num_items, items);
               if (!cluster->distances)
                    goto cleanup;
               cluster->num_items = num_items;
               cluster->num_clusters = 0;
               if (add_leaves(cluster, items))
                    merge_clusters(cluster);
               else
                    goto cleanup;
          } else {
               fprintf(stderr, "Failed to allocate cluster nodes.\n");
               goto cleanup;
          }
     } else
          fprintf(stderr, "Failed to allocate cluster.\n");
     goto done;

cleanup:
     free_cluster(cluster);
     cluster = NULL;

done:
     return cluster;
}

void print_cluster_node(cluster_node_t *node) {
     proximity_node_t *t;
     if (node->label)
          fprintf(stdout, "%s\n", node->label);
     t = node->proximity;
     while (t) {
          fprintf(stdout, "(%d, %f), ", t->target, t->distance);
          t = t->next;
     }
     fprintf(stdout, "\n");
}

void print_cluster(cluster_t *cluster) {
     int i;
     for (i = 0; i < cluster->num_clusters; ++i)
          print_cluster_node(&cluster->nodes[i]);
}

int main(int argc, char **argv) {
     item_t items[] = {
          { .label = "A", .x = 0.5, .y = 0.5 },
          { .label = "B", .x = 5.5, .y = 0.5 },
          { .label = "C", .x = 5.5, .y = 5.5 },
          { .label = "D", .x = 0.5, .y = 5.5 }
     };
     cluster_t *cluster = agglomerate(sizeof(items) / sizeof(items[0]), items);
     print_cluster(cluster);
     free_cluster(cluster);
     return 0;
}
