#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "utils.h"
#define _FILE_OFFSET_BITS 64

#include "output_datatype.h"

#ifndef ADD_DIFF_TIME
#define ADD_DIFF_TIME(tbegin, tend)  ((tend.tv_sec - tstart.tv_sec) + 1e-6*(tend.tv_usec - tstart.tv_usec))
#endif

extern int CTREES_UPID_FEATURE;

/* Make sure LOCATIONS_FILENAME is a multiple of 8 */
#define LOCATIONS_FILENAME_SIZE  (32)
struct locations{
    int64_t forestid; /* forest id for the tree */
    int64_t tree_root; /* tree root id (id for z=0 halo) */
    int64_t fileid;/* fileid, 1d index for the 3-d decomposition by consistent-trees */
    int64_t offset;/* the offset in tree_* file. */
    int64_t bytes;/* number of bytes in the file. MUST BE SIGNED type */
    char filename[LOCATIONS_FILENAME_SIZE];/* filename where the tree is written */
};

struct forest_info{
    int64_t nforests;
    int64_t *forestid;
    int64_t *fileid;
    int64_t *num_trees;
    int64_t *num_halos;
    int64_t *num_ascii_bytes;
    int64_t *num_binary_bytes;
    int64_t *offset;
    char (*filename)[LOCATIONS_FILENAME_SIZE];
};    

struct additional_info{
    int64_t id;
    int64_t pid;
    int64_t upid;
    double desc_scale;
    int64_t descid;
    double scale;
};


int64_t read_forests(const char *filename, const char *output_dir, int64_t **f, int64_t **t);
int64_t read_locations(const char *filename, const int64_t ntrees, struct locations *l, int64_t *num_files, int64_t *box_div);
void sort_locations_file_offset(const int64_t ntrees, struct locations *locations);
int run_checks_on_new_locations(const struct locations *new_locations, const struct locations *locations, const int64_t ntrees);
void sort_forests(const int64_t ntrees, int64_t *forests, int64_t *tree_roots);
void assign_forest_ids(const int64_t ntrees, struct locations *locations, int64_t *forests, int64_t *tree_roots);
struct forest_info * assign_trees_in_forest_to_same_file(const int64_t ntrees, struct locations *locations, struct locations *output_locations, const int64_t nfiles);
int64_t read_tree_into_forest(int64_t *nhalos_allocated, struct output_dtype **source_forest, const int64_t forest_offset,
                              struct additional_info **source_info, 
#ifdef USE_FGETS                              
                              FILE *fp,
                              int64_t offset, 
#else
                              int fd,
                              off_t offset, 
#endif
                              const size_t bytes,
                              const float inv_part_mass);


int64_t compute_numbytes_with_off(const int64_t off, const int64_t start);
int64_t compute_numbytes(FILE *fp, const int64_t start);

int64_t write_forests_and_locations(const char *filename, const int64_t ntrees, const struct locations *locations);
int64_t find_fof_halo(const int64_t totnhalos, const struct additional_info *info, const int start_loc, const int64_t upid, int verbose, int64_t calldepth);
int fix_upid(const int64_t totnhalos, struct output_dtype *forest, struct additional_info *info, int *interrupted, const int verbose);
void assign_mergertree_indices(const int64_t totnhalos, struct output_dtype *forest, struct additional_info *info, const int max_snapnum);
int64_t ** calculate_forest_info_per_file(const int64_t nfiles, int64_t *totnforests_per_file, const struct forest_info *forest_info);
int fix_flybys(const int64_t totnhalos, struct output_dtype *forest, struct additional_info *info, int verbose);
void validate_fields(const int64_t totnhalos, const struct output_dtype *forest, const struct additional_info *info, const int max_snapnum);
    
    
