/* Copyright 2014 Gagarine Yaikhom (MIT License) */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NOT_USED  0 /* node is currently not used */
#define LEAF_NODE 1 /* node contains a leaf node */
#define A_MERGER  2 /* node contains a merged pair of root clusters */
#define MAX_LABEL_LEN 16

#define SINGLE_LINKAGE   's' /* choose minimum distance */
#define COMPLETE_LINKAGE 'c' /* choose maximum distance */
#define AVERAGE_LINKAGE  'a' /* choose average distance */

typedef struct cluster_s cluster_t;
typedef struct cluster_node_s cluster_node_t;
typedef struct neighbour_s neighbour_t;
typedef struct item_s item_t;

float (*distance_fptr)(float **, const int *, const int *, int, int);

struct cluster_s {
     int num_items; /* number of items that was clustered */
     int num_clusters; /* current number of root clusters */
     int num_nodes; /* number of leaf and merged clusters */
     cluster_node_t *nodes; /* leaf and merged clusters */
     float **distances; /* distance between leaves */
};

struct cluster_node_s {
     int type; /* type of the cluster node */
     int is_root; /* true if cluster hasn't merged with another */
     int height; /* height of node from the bottom */
     char *label; /* label of a leaf node */
     int *merged; /* indexes of root clusters merged */
     int num_items; /* number of leaf nodes inside new cluster */
     int *items; /* array of leaf nodes indices inside merged clusters */
     neighbour_t *neighbours; /* sorted linked list of distances to roots */
};

struct neighbour_s {
     int target; /* the index of cluster node representing neighbour */
     float distance; /* distance between the nodes */
     neighbour_t *next, *prev; /* linked list entries */
};

struct item_s {
     float x, y; /* (x, y) coordinate of the input data point */
     char label[MAX_LABEL_LEN]; /* label of the input data point */
};

cluster_node_t *add_leaf(cluster_t *cluster, const char *label) {
     cluster_node_t *t = &(cluster->nodes[cluster->num_nodes]);
     t->label = (char *) calloc(strlen(label), sizeof(char));
     if (t->label) {
          t->items = (int *) calloc(1, sizeof(int));
          if (t->items) {
               strlcpy(t->label, label, sizeof(label));
               t->type = LEAF_NODE;
               t->is_root = 1;
               t->height = 0;
               t->num_items = 1;
               t->items[0] = cluster->num_nodes++;
               cluster->num_clusters++;
          } else {
               fprintf(stderr, "Failed to allocate items.\n");
               free(t->label);
               t = NULL;
          }
     } else {
          fprintf(stderr, "Failed to allocate label.\n");
          t = NULL;
     }
     return t;
}

void free_neighbours(neighbour_t *node) {
     neighbour_t *t;
     while (node) {
          t = node->next;
          free(node);
          node = t;
     }
}

void free_cluster(cluster_t * cluster) {
     if (cluster) {
          if (cluster->nodes) {
               for (int i = 0; i < cluster->num_nodes; ++i) {
                    cluster_node_t node = cluster->nodes[i];
                    if (node.label)
                         free(node.label);
                    if (node.neighbours)
                         free_neighbours(node.neighbours);
               }
               free(cluster->nodes);
          }
          if (cluster->distances)
               free(cluster->distances);
          free(cluster);
     }
}

void insert_before(neighbour_t *current, neighbour_t *neighbours,
                   cluster_node_t *node) {
     neighbours->next = current;
     if (current->prev) {
          current->prev->next = neighbours;
          neighbours->prev = current->prev;
     } else
          node->neighbours = neighbours;
     current->prev = neighbours;
}

void insert_after(neighbour_t *current, neighbour_t *neighbours) {
     neighbours->prev = current;
     current->next = neighbours;
}

void insert_sorted(cluster_node_t *node, neighbour_t *neighbours) {
     neighbour_t *temp = node->neighbours;
     while (temp->next) {
          if (temp->distance >= neighbours->distance) {
               insert_before(temp, neighbours, node);
               return;
          }
          temp = temp->next;
     }
     if (neighbours->distance < temp->distance)
          insert_before(temp, neighbours, node);
     else
          insert_after(temp, neighbours);
}

/**
 * Find pair-wise shortest distance between cluster leaves.
 */
float single_linkage_distance(float **distances, const int a[],
                              const int b[], int m, int n) {
     float min = FLT_MAX, d;
     for (int i = 0; i < m; ++i)
          for (int j = 0; j < n; ++j) {
               d = distances[a[i]][b[j]];
               if (d < min)
                    min = d;
          }
     return min;
}

/**
 * Find pair-wise longest distance between cluster leaves.
 */
float complete_linkage_distance(float **distances, const int a[],
                                const int b[], int m, int n) {
     float d, max = 0.0 /* assuming distances are positive */;
     for (int i = 0; i < m; ++i)
          for (int j = 0; j < n; ++j) {
               d = distances[a[i]][b[j]];
               if (d > max)
                    max = d;
          }
     return max;
}

/**
 * Find average of pair-wise distances between cluster leaves.
 */
float average_linkage_distance(float **distances, const int a[],
                               const int b[], int m, int n) {
     float total = 0.0;
     for (int i = 0; i < m; ++i)
          for (int j = 0; j < n; ++j)
               total += distances[a[i]][b[j]];
     return total / (m * n);
}

float get_distance(cluster_t *cluster, int index, int target) {
     /* if both are leaves, just use the distances matrix */
     if (index < cluster->num_items && target < cluster->num_items)
          return cluster->distances[index][target];
     else
          return distance_fptr(cluster->distances,
                               cluster->nodes[index].items,
                               cluster->nodes[target].items,
                               cluster->nodes[index].num_items,
                               cluster->nodes[target].num_items);
}

neighbour_t *add_neighbour(cluster_t *cluster, int index, int target) {
     neighbour_t *neighbour = (neighbour_t *)
          calloc(1, sizeof(neighbour_t));
     if (neighbour) {
          neighbour->target = target;
          neighbour->distance = get_distance(cluster, index, target);
          cluster_node_t *node = &(cluster->nodes[index]);
          if (node->neighbours)
               insert_sorted(node, neighbour);
          else
               node->neighbours = neighbour;
     } else
          fprintf(stderr, "Failed to allocate neighbour node.\n");
     return neighbour;
}

cluster_t *update_neighbours(cluster_t *cluster, int index) {
     cluster_node_t node = cluster->nodes[index];
     if (node.type == NOT_USED) {
          fprintf(stderr, "Supplied node with index %d is invalid\n", index);
          cluster = NULL;
     } else {
          int root_clusters_seen = 1, target = index;
          while (root_clusters_seen < cluster->num_clusters) {
               cluster_node_t temp = cluster->nodes[--target];
               if (temp.type == NOT_USED) {
                    fprintf(stderr, "Empty node found at %d.\n", index);
                    cluster = NULL;
                    break;
               }
               if (temp.is_root) {
                    ++root_clusters_seen;
                    add_neighbour(cluster, index, target);
               }
          }
     }
     return cluster;
}

float euclidean_distance(item_t a, item_t b) {
     return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
}

void fill_euclidean_distances(float **matrix, int num_items,
                             const item_t items[]) {
     for (int i = 0; i < num_items; ++i)
          for (int j = 0; j < num_items; ++j) {
               matrix[i][j] = euclidean_distance(items[i], items[j]);
               matrix[j][i] = matrix[i][j];
          }
}

float **generate_distance_matrix(int num_items, const item_t items[]) {
     float **matrix = (float **) calloc(num_items, sizeof(float *));
     if (matrix) {
          for (int i = 0; i < num_items; ++i) {
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
               fill_euclidean_distances(matrix, num_items, items);
     } else
          fprintf(stderr, "Failed to allocate distance matrix.\n");
     return matrix;
}

cluster_t *add_leaves(cluster_t *cluster, item_t *items) {
     for (int i = 0; i < cluster->num_items; ++i) {
          if (add_leaf(cluster, items[i].label))
               update_neighbours(cluster, i);
          else {
               cluster = NULL;
               break;
          }
     }
     free(items);
     return cluster;
}

cluster_node_t *merge(cluster_t *cluster, int first_idx) {
     int new_idx = cluster->num_nodes;
     cluster_node_t *node = &(cluster->nodes[new_idx]);
     node->merged = (int *) calloc(2, sizeof(int));
     if (node->merged) {
          cluster_node_t *first = &(cluster->nodes[first_idx]);
          int second_idx = first->neighbours->target;
          cluster_node_t *second = &(cluster->nodes[second_idx]);
        
          /* expand and combine leaf items that was merged */
          node->num_items = first->num_items + second->num_items;
          node->items = (int *) calloc(node->num_items, sizeof(int));
          if (node->items) {
               /* what was merged? */
               node->type = A_MERGER;
               node->merged[0] = first_idx;
               node->merged[1] = second_idx;
            
               /* copy leaf indexes from merged clusters */
               int j = 0;
               for (int i = 0; i < first->num_items; ++i)
                    node->items[j++] = first->items[i];
               for (int i = 0; i < second->num_items; ++i)
                    node->items[j++] = second->items[i];
            
               /* update root clusters */
               node->is_root = 1;
               first->is_root = 0;
               second->is_root = 0;
            
               /* set node height of new cluster */
               node->height = first->height;
               if (node->height < second->height)
                    node->height = second->height;
               node->height++;
            
               cluster->num_nodes++;
               cluster->num_clusters--;
               update_neighbours(cluster, new_idx);
          }
     } else {
          fprintf(stderr, "Failed to allocate merge array.\n");
          node = NULL;
     }
     return node;
}

cluster_t *merge_clusters(cluster_t *cluster) {
     float best_distance = 0.0;
     while (cluster->num_clusters > 1) {
          int root_clusters_seen = 0;
          int j = cluster->num_nodes - 1; /* traverse hierarchy top-down */
          int best_candidate = -1;
          while (root_clusters_seen < cluster->num_clusters) {
               cluster_node_t node = cluster->nodes[j];
               if (node.type != NOT_USED && node.is_root) {
                    ++root_clusters_seen;
                    if (node.neighbours)
                         if (best_candidate == -1 ||
                             node.neighbours->distance < best_distance) {
                              best_candidate = j;
                              best_distance = node.neighbours->distance;
                         }
               }
               --j;
          }
          merge(cluster, best_candidate);
     }
     return cluster;
}

cluster_t *agglomerate(int num_items, item_t *items) {
     cluster_t *cluster = (cluster_t *) calloc(1, sizeof(cluster_t));
     if (cluster) {
          cluster->nodes = (cluster_node_t *)
               calloc(2 * num_items - 1, sizeof(cluster_node_t));
          if (cluster->nodes) {
               cluster->distances = generate_distance_matrix(num_items, items);
               if (!cluster->distances)
                    goto cleanup;
               cluster->num_items = num_items;
               cluster->num_nodes = 0;
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

void print_cluster_node(cluster_t *cluster, int index) {
     neighbour_t *t;
     cluster_node_t node = cluster->nodes[index];
     fprintf(stdout, "Node %d (height: %d)\n", index, node.height);
     if (node.label)
          fprintf(stdout, "\tLeaf: %s\n", node.label);
     else
          fprintf(stdout, "\tMerged: %d, %d\n", node.merged[0], node.merged[1]);
     fprintf(stdout, "\tItems: ");
     if (node.num_items > 0) {
          fprintf(stdout, "%s", cluster->nodes[node.items[0]].label);
          for (int i = 1; i < node.num_items; ++i)
               fprintf(stdout, ", %s", cluster->nodes[node.items[i]].label);
     }
     fprintf(stdout, "\n\tNeighbours: ");
     t = node.neighbours;
     while (t) {
          fprintf(stdout, "\n\t\t%2d: %5.3f", t->target, t->distance);
          t = t->next;
     }
     fprintf(stdout, "\n");
}

void print_cluster(cluster_t *cluster) {
     int i;
     for (i = 0; i < cluster->num_nodes; ++i)
          print_cluster_node(cluster, i);
}

int process_input(item_t **items, const char *fname) {
     int count = 0;
     FILE *f = fopen(fname, "r");
     if (f) {
          char linkage_type;
          fscanf(f, "%c\n", &linkage_type);
          switch (linkage_type) {
          case AVERAGE_LINKAGE:
               distance_fptr = average_linkage_distance;
               break;
          case COMPLETE_LINKAGE:
               distance_fptr = complete_linkage_distance;
               break;
          case SINGLE_LINKAGE:
          default: distance_fptr = single_linkage_distance;
          }
        
          fscanf(f, "%d\n", &count);
          if (count) {
               item_t *temp = (item_t *) calloc(count, sizeof(item_t));
               if (temp) {
                    for (int i = 0; i < count; ++i)
                         if (!fscanf(f, "%[^|]| %f %f\n", temp[i].label,
                                     &(temp[i].x), &(temp[i].y))) {
                              fprintf(stderr, "Failed to read input.\n");
                              break;
                         }
                    *items = temp;
               } else
                    fprintf(stderr, "Failed to allocate items array.\n");
          }
     } else
          fprintf(stderr, "Failed to open input file.\n");
     return count;
}

int main(int argc, char **argv) {
     if (argc != 2) {
          fprintf(stderr, "Usage: %s <input file>\n", argv[0]);
          exit(1);
     } else {
          item_t *items = NULL;
          int num_items = process_input(&items, argv[1]);
          if (num_items) {
               cluster_t *cluster = agglomerate(num_items, items);
               print_cluster(cluster);
               free_cluster(cluster);
          }
     }
     return 0;
}
