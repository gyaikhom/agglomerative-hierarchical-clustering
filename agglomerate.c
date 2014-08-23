/**
 * Copyright 2014 Gagarine Yaikhom (MIT License)
 *
 * Implements Agglomerative Hierarchical Clustering algorithm.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NOT_USED  0 /* node is currently not used */
#define LEAF_NODE 1 /* node contains a leaf node */
#define A_MERGER  2 /* node contains a merged pair of root clusters */
#define MAX_LABEL_LEN 16

#define AVERAGE_LINKAGE  'a' /* choose average distance */
#define CENTROID_LINKAGE 't' /* choose distance between cluster centroids */
#define COMPLETE_LINKAGE 'c' /* choose maximum distance */
#define SINGLE_LINKAGE   's' /* choose minimum distance */

#define alloc_mem(N, T) (T *) calloc(N, sizeof(T))
#define alloc_fail(M) fprintf(stderr,                                   \
                              "Failed to allocate memory for %s.\n", M)
#define read_fail(M) fprintf(stderr, "Failed to read %s from file.\n", M)
#define invalid_node(I) fprintf(stderr,                                 \
                                "Invalid cluster node at index %d.\n", I)

typedef struct cluster_s cluster_t;
typedef struct cluster_node_s cluster_node_t;
typedef struct neighbour_s neighbour_t;
typedef struct item_s item_t;

float (*distance_fptr)(float **, const int *, const int *, int, int);

typedef struct coord_s {
        float x, y;
} coord_t;

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
        coord_t centroid; /* centroid of this cluster */
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
        coord_t coord; /* coordinate of the input data point */
        char label[MAX_LABEL_LEN]; /* label of the input data point */
};

float euclidean_distance(const coord_t *a, const coord_t *b) {
        return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2));
}

cluster_node_t *add_leaf(cluster_t *cluster, const item_t *item) {
        cluster_node_t *t = &(cluster->nodes[cluster->num_nodes]);
        int len = strlen(item->label) + 1;
        t->label = alloc_mem(len, char);
        if (t->label) {
                t->items = alloc_mem(1, int);
                if (t->items) {
                        strncpy(t->label, item->label, len);
                        t->centroid = item->coord;
                        t->type = LEAF_NODE;
                        t->is_root = 1;
                        t->height = 0;
                        t->num_items = 1;
                        t->items[0] = cluster->num_nodes++;
                        cluster->num_clusters++;
                } else {
                        alloc_fail("node items");
                        free(t->label);
                        t = NULL;
                }
        } else {
                alloc_fail("node label");
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
                                if (node.merged)
                                        free(node.merged);
                                if (node.items)
                                        free(node.items);
                                if (node.neighbours)
                                        free_neighbours(node.neighbours);
                        }
                        free(cluster->nodes);
                }
                if (cluster->distances) {
                        for (int i = 0; i < cluster->num_items; ++i)
                                free(cluster->distances[i]);
                        free(cluster->distances);
                }
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

float single_linkage(float **distances, const int a[],
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

float complete_linkage(float **distances, const int a[],
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

float average_linkage(float **distances, const int a[],
                      const int b[], int m, int n) {
        float total = 0.0;
        for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                        total += distances[a[i]][b[j]];
        return total / (m * n);
}

float centroid_linkage(float **distances, const int a[],
                       const int b[], int m, int n) {
        return 0; /* empty function */
}

float get_distance(cluster_t *cluster, int index, int target) {
        /* if both are leaves, just use the distances matrix */
        if (index < cluster->num_items && target < cluster->num_items)
                return cluster->distances[index][target];
        else {
                cluster_node_t *a = &(cluster->nodes[index]);
                cluster_node_t *b = &(cluster->nodes[target]);
                if (distance_fptr == centroid_linkage)
                        return euclidean_distance(&(a->centroid),
                                                  &(b->centroid));
                else return distance_fptr(cluster->distances,
                                          a->items, b->items,
                                          a->num_items, b->num_items);
        }
}

neighbour_t *add_neighbour(cluster_t *cluster, int index, int target) {
        neighbour_t *neighbour = alloc_mem(1, neighbour_t);
        if (neighbour) {
                neighbour->target = target;
                neighbour->distance = get_distance(cluster, index, target);
                cluster_node_t *node = &(cluster->nodes[index]);
                if (node->neighbours)
                        insert_sorted(node, neighbour);
                else
                        node->neighbours = neighbour;
        } else
                alloc_fail("neighbour node");
        return neighbour;
}

cluster_t *update_neighbours(cluster_t *cluster, int index) {
        cluster_node_t node = cluster->nodes[index];
        if (node.type == NOT_USED) {
                invalid_node(index);
                cluster = NULL;
        } else {
                int root_clusters_seen = 1, target = index;
                while (root_clusters_seen < cluster->num_clusters) {
                        cluster_node_t temp = cluster->nodes[--target];
                        if (temp.type == NOT_USED) {
                                invalid_node(index);
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

void fill_euclidean_distances(float **matrix, int num_items,
                              const item_t items[]) {
        for (int i = 0; i < num_items; ++i)
                for (int j = 0; j < num_items; ++j) {
                        matrix[i][j] =
                                euclidean_distance(&(items[i].coord),
                                                   &(items[j].coord));
                        matrix[j][i] = matrix[i][j];
                }
}

float **generate_distance_matrix(int num_items, const item_t items[]) {
        float **matrix = alloc_mem(num_items, float *);
        if (matrix) {
                for (int i = 0; i < num_items; ++i) {
                        matrix[i] = alloc_mem(num_items, float);
                        if (!matrix[i]) {
                                alloc_fail("distance matrix row");
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
                alloc_fail("distance matrix");
        return matrix;
}

cluster_t *add_leaves(cluster_t *cluster, item_t *items) {
        for (int i = 0; i < cluster->num_items; ++i) {
                if (add_leaf(cluster, &items[i]))
                        update_neighbours(cluster, i);
                else {
                        cluster = NULL;
                        break;
                }
        }
        return cluster;
}

void print_cluster_items(cluster_t *cluster, int index) {
        cluster_node_t *node = &(cluster->nodes[index]);
        fprintf(stdout, "Items: ");
        if (node->num_items > 0) {
                fprintf(stdout, "%s", cluster->nodes[node->items[0]].label);
                for (int i = 1; i < node->num_items; ++i)
                        fprintf(stdout, ", %s",
                                cluster->nodes[node->items[i]].label);
        }
        fprintf(stdout, "\n");
}

void print_cluster_node(cluster_t *cluster, int index) {
        cluster_node_t *node = &(cluster->nodes[index]);
        fprintf(stdout, "Node %d - height: %d, centroid: (%5.3f, %5.3f)\n",
                index, node->height, node->centroid.x, node->centroid.y);
        if (node->label)
                fprintf(stdout, "\tLeaf: %s\n\t", node->label);
        else
                fprintf(stdout, "\tMerged: %d, %d\n\t",
                        node->merged[0], node->merged[1]);
        print_cluster_items(cluster, index);
        fprintf(stdout, "\tNeighbours: ");
        neighbour_t *t = node->neighbours;
        while (t) {
                fprintf(stdout, "\n\t\t%2d: %5.3f", t->target, t->distance);
                t = t->next;
        }
        fprintf(stdout, "\n");
}

cluster_node_t *merge(cluster_t *cluster, int first_idx, int second_idx) {
        int new_idx = cluster->num_nodes;
        cluster_node_t *node = &(cluster->nodes[new_idx]);
        node->merged = alloc_mem(2, int);
        if (node->merged) {
                cluster_node_t *to_merge[2] = {
                        &(cluster->nodes[first_idx]),
                        &(cluster->nodes[second_idx])
                };
        
                /* expand and combine leaf items that was merged */
                node->num_items = to_merge[0]->num_items + to_merge[1]->num_items;
                node->items = alloc_mem(node->num_items, int);
                if (node->items) {
                        /* what was merged? */
                        node->type = A_MERGER;
                        node->is_root = 1;
                        node->merged[0] = first_idx;
                        node->merged[1] = second_idx;
                        node->height = -1;

                        /* copy leaf indexes from merged clusters */
                        int k = 0, idx;
                        coord_t centroid = { .x = 0.0, .y = 0.0 };
                        for (int i = 0; i < 2; ++i) {
                                cluster_node_t *t = to_merge[i];
                                t->is_root = 0; /* no longer root after merger */
                                if (node->height == -1 || node->height < t->height) 
                                        node->height = t->height;
                                for (int j = 0; j < t->num_items; ++j) {
                                        idx = t->items[j];
                                        node->items[k++] = idx;

                                        /* for calculating centroid of new cluster */
                                        centroid.x += cluster->nodes[idx].centroid.x;
                                        centroid.y += cluster->nodes[idx].centroid.y;
                                }
                        }
                        /* calculate centroid */
                        node->centroid.x = centroid.x / k;
                        node->centroid.y = centroid.y / k;

                        node->height++;
                        cluster->num_nodes++;
                        cluster->num_clusters--;
                        update_neighbours(cluster, new_idx);
                } else {
                        alloc_fail("array of merged items");
                        free(node->merged);
                        node = NULL;
                }
        } else {
                alloc_fail("array of merged nodes");
                node = NULL;
        }
        return node;
}

int find_clusters_to_merge(cluster_t *cluster, int *first_idx, int *second_idx) {
        float best_distance = 0.0;
        int root_clusters_seen = 0;
        int j = cluster->num_nodes; /* traverse hierarchy top-down */
        *first_idx = -1;
        while (root_clusters_seen < cluster->num_clusters) {
                cluster_node_t node = cluster->nodes[--j];
                if (node.type == NOT_USED || !node.is_root)
                        continue;
                ++root_clusters_seen;
                neighbour_t *t = node.neighbours;
                while (t) {
                        if (cluster->nodes[t->target].is_root) {
                                if (*first_idx == -1 ||
                                    t->distance < best_distance) {
                                        *first_idx = j;
                                        *second_idx = t->target;
                                        best_distance = t->distance;
                                }
                                break;
                        }
                        t = t->next;
                }
        }
        return *first_idx;
}

cluster_t *merge_clusters(cluster_t *cluster) {
        int first, second;
        while (cluster->num_clusters > 1) {
                if (find_clusters_to_merge(cluster, &first, &second) != -1)
                        merge(cluster, first, second);
        }
        return cluster;
}

cluster_t *agglomerate(int num_items, item_t *items) {
        cluster_t *cluster = alloc_mem(1, cluster_t);
        if (cluster) {
                cluster->nodes = alloc_mem(2 * num_items - 1, cluster_node_t);
                if (cluster->nodes) {
                        cluster->distances =
                                generate_distance_matrix(num_items, items);
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
                        alloc_fail("cluster nodes");
                        goto cleanup;
                }
        } else
                alloc_fail("cluster");
        goto done;
    
cleanup:
        free_cluster(cluster);
        cluster = NULL;
    
done:
        return cluster;
}

void print_cluster(cluster_t *cluster) {
        for (int i = 0; i < cluster->num_nodes; ++i)
                print_cluster_node(cluster, i);
}

int read_items_from_file(item_t **items, FILE *f) {
        int count, r;
        r = fscanf(f, "%5d\n", &count);
        if (r == 0) {
                read_fail("number of lines");
                return 0;
        }
        if (count) {
                item_t *t = alloc_mem(count, item_t);
                if (t) {
                        for (int i = 0; i < count; ++i) {
                                r = fscanf(f, "%[^|]| %10f %10f\n",
                                           t[i].label, &(t[i].coord.x), &(t[i].coord.y));
                                if (r == 0) {
                                        read_fail("item line");
                                        free(t);
                                        return 0;
                                }
                        }
                        *items = t;
                } else
                        alloc_fail("items array");
        }
        return count;
}

int process_input(item_t **items, const char *fname) {
        int count = 0;
        FILE *f = fopen(fname, "r");
        if (f) {
                count = read_items_from_file(items, f);
                fclose(f);
        } else
                fprintf(stderr, "Failed to open input file %s.\n", fname);
        return count;
}

void get_k_clusters(cluster_t *cluster, int k) {
        if (k < 1)
                return;
        if (k > cluster->num_items)
                k = cluster->num_items;

        int i = cluster->num_nodes - 1;
        int nodes_to_discard = cluster->num_nodes - k + 1;
        while (k) {
                if (i < nodes_to_discard) {
                        print_cluster_items(cluster, i);
                        --k;
                } else {
                        cluster_node_t *node = &(cluster->nodes[i]);
                        if (node->type == A_MERGER) {
                                for (int j = 0; j < 2; ++j) {
                                        int t = node->merged[j];
                                        if (t < nodes_to_discard) {
                                                print_cluster_items(cluster, t);
                                                --k;
                                        }
                                }
                        }
                }
                --i;
        }
}

void set_linkage(char linkage_type) {
        switch (linkage_type) {
        case AVERAGE_LINKAGE:
                distance_fptr = average_linkage;
                break;
        case COMPLETE_LINKAGE:
                distance_fptr = complete_linkage;
                break;
        case CENTROID_LINKAGE:
                distance_fptr = centroid_linkage;
                break;
        case SINGLE_LINKAGE:
        default: distance_fptr = single_linkage;
        }
}

int main(int argc, char **argv) {
        if (argc != 4) {
                fprintf(stderr, "Usage: %s <input file> <num clusters> "
                        "<linkage type>\n", argv[0]);
                exit(1);
        } else {
                item_t *items = NULL;
                int num_items = process_input(&items, argv[1]);
                set_linkage(argv[3][0]);
                if (num_items) {
                        cluster_t *cluster = agglomerate(num_items, items);
                        free(items);

                        fprintf(stdout, "CLUSTER HIERARCHY\n--------------------\n");
                        print_cluster(cluster);
            
                        int k = atoi(argv[2]);
                        fprintf(stdout, "\n\n%d CLUSTERS\n--------------------\n", k);
                        get_k_clusters(cluster, k);
                        free_cluster(cluster);
                }
        }
        return 0;
}
