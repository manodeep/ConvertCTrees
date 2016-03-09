#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>

#define _FILE_OFFSET_BITS 64

#include "utils.h"
#include "sglib.h"
#include "progressbar.h"

#include "output_datatype.h"

#ifndef ADD_DIFF_TIME
#define ADD_DIFF_TIME(tbegin, tend)  ((tend.tv_sec - tstart.tv_sec) + 1e-6*(tend.tv_usec - tstart.tv_usec))
#endif

/* #define FIELD_CHECKER(field, Nmax, fieldname)   (XASSERT((field == -1) || (field >=0 && field < Nmax), fieldname " = %"PRId64" must be " */

static int CTREES_UPID_BUG = 0;

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
    int64_t *num_bytes;
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


void usage(int argc, char **argv)
{
    (void) argc;
    fprintf(stderr,"USAGE: %s <input consistent-trees directory> <output LHALOTREE directory> <particle mass (10^10 Msun/h units) >\n",
            argv[0]);
}    

int64_t read_forests(const char *filename, int64_t **f, int64_t **t)
{
    char buffer[MAXBUFSIZE];
    const char comment = '#';
    /* By passing the comment character, getnumlines
       will return the actual number of lines, ignoring
       the first header line. 
     */

    const int64_t ntrees = getnumlines(filename, comment);
    *f = my_malloc(sizeof(int64_t), ntrees);
    *t = my_malloc(sizeof(int64_t), ntrees);

    int64_t *forests    = *f;
    int64_t *tree_roots = *t;
    
    int64_t ntrees_found = 0;
    FILE *fp = my_fopen(filename, "r");
    while(fgets(buffer, MAXBUFSIZE, fp) != NULL) {
        if(buffer[0] == comment) {
            continue;
        } else {
            const int nitems_expected = 2;
            XASSERT(ntrees_found < ntrees,
                    "ntrees=%"PRId64" should be less than ntrees_found=%"PRId64"\n", ntrees, ntrees_found);            
            int nitems = sscanf(buffer, "%"SCNd64" %"SCNd64, &(tree_roots[ntrees_found]), &(forests[ntrees_found]));
            XASSERT(nitems == nitems_expected,
                    "Expected to parse %d long integers but found `%s' in the buffer. nitems = %d \n",
                    nitems_expected, buffer, nitems);
            ntrees_found++;
        }
    }
    XASSERT(ntrees == ntrees_found,
            "ntrees=%"PRId64" should be equal to ntrees_found=%"PRId64"\n", ntrees, ntrees_found);
    fclose(fp);

    return ntrees_found;
}    


int64_t read_locations(const char *filename, const int64_t ntrees, struct locations *l, int *num_files, int *box_div)
{
    char buffer[MAXBUFSIZE];
    int max_fileid = 0;
    const char comment = '#';
    /* By passing the comment character, getnumlines
       will return the actual number of lines, ignoring
       the first header line. 
     */

    struct locations *locations = l;

    int64_t ntrees_found = 0;
    FILE *fp = my_fopen(filename, "r");
    while(fgets(buffer, MAXBUFSIZE, fp) != NULL) {
        if(buffer[0] == comment) {
            continue;
        } else {
            const int nitems_expected = 4;
            char linebuf[MAXLEN];
            XASSERT(ntrees_found < ntrees,
                    "ntrees=%"PRId64" should be less than ntrees_found=%"PRId64"\n",
                    ntrees, ntrees_found);            
            int nitems = sscanf(buffer, "%"SCNd64" %"SCNd64 " %"SCNd64 "%s", &(locations[ntrees_found].tree_root),
                                &(locations[ntrees_found].fileid), &(locations[ntrees_found].offset), linebuf);

            XASSERT(locations[ntrees_found].offset >= 0,
                    "offset=%"PRId64" for ntree =%"PRId64" must be positive.\nFile = `%s'\nbuffer = `%s'\n",
                    locations[ntrees_found].offset,ntrees_found,filename, buffer);

            
            /* The filename is separated out to save memory but I want to ensure that the actual filename does
               not get truncated. The filename field might actually be removed later. */
            my_snprintf(locations[ntrees_found].filename, LOCATIONS_FILENAME_SIZE, "%s", linebuf);
            XASSERT(nitems == nitems_expected,
                    "Expected to parse two long integers but found `%s' in the buffer\n",
                    buffer);
            ntrees_found++;
        }
    }
    XASSERT(ntrees == ntrees_found, "ntrees=%"PRId64" should be equal to ntrees_found=%"PRId64"\n", ntrees, ntrees_found);
    fclose(fp);

    for(int i=0;i<ntrees_found;i++){
        if (locations[i].fileid > max_fileid) {
            max_fileid = locations[i].fileid;
        }
    }

    /* number of files is one greater from 0 based indexing of C files */
    *num_files = max_fileid + 1;
    const int box_divisions = (int) round(cbrt(*num_files));
    const int box_cube = box_divisions * box_divisions * box_divisions;
    XASSERT( (box_cube) == (*num_files),
             "box_divisions^3=%d should be equal to nfiles=%d\n",
             box_cube, *num_files);
    *box_div = box_divisions;
    return ntrees_found;
}    


void sort_forests(const int64_t ntrees, int64_t *forests, int64_t *tree_roots)
{

#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64_t, tree_roots,i,j); \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64_t, forests, i, j) }
    
    SGLIB_ARRAY_QUICK_SORT(int64_t, tree_roots, ntrees, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER

}    

int compare_locations_tree_roots(const void *l1, const void *l2)
{
    const struct locations *aa = (const struct locations *) l1;
    const struct locations *bb = (const struct locations *) l2;
    return (aa->tree_root < bb->tree_root) ? -1:1;
}

int compare_locations_fid(const void *l1, const void *l2)
{
    const struct locations *aa = (const struct locations *) l1;
    const struct locations *bb = (const struct locations *) l2;
    return (aa->forestid < bb->forestid) ? -1:1;
}

int compare_locations_file_offset(const void *l1, const void *l2)
{
    const struct locations *aa = (const struct locations *) l1;
    const struct locations *bb = (const struct locations *) l2;

    const int file_id_cmp = (aa->fileid == bb->fileid) ? 0:((aa->fileid < bb->fileid) ? -1:1);
    if(file_id_cmp == 0) {
        /* trees are in same file -> sort by offset */
        return (aa->offset < bb->offset) ? -1:1;
    } else {
        return file_id_cmp;
    }
        
    return 0;
}

int compare_locations_fid_file_offset(const void *l1, const void *l2)
{
    const struct locations *aa = (const struct locations *) l1;
    const struct locations *bb = (const struct locations *) l2;

    if(aa->forestid != bb->forestid) {
        return (aa->forestid < bb->forestid) ? -1:1;
    } else {
        /* The trees are in the same forest. Check filename */
        const int file_id_cmp = (aa->fileid == bb->fileid) ? 0:((aa->fileid < bb->fileid) ? -1:1);
        if(file_id_cmp == 0) {
            /* trees are in same file -> sort by offset */
            return (aa->offset < bb->offset) ? -1:1;
        } else {
            return file_id_cmp;
        }
    }
        
    return 0;
}

void sort_locations_on_treeroot(const int64_t ntrees, struct locations *locations)
{
    qsort(locations, ntrees, sizeof(*locations), compare_locations_tree_roots);
}    

void sort_locations_file_offset(const int64_t ntrees, struct locations *locations)
{
    qsort(locations, ntrees, sizeof(*locations), compare_locations_file_offset);
}

void sort_locations_on_fid(const int64_t ntrees, struct locations *locations)
{
    qsort(locations, ntrees, sizeof(*locations), compare_locations_fid);
}    

void sort_locations_on_fid_file_offset(const int64_t ntrees, struct locations *locations)
{
    qsort(locations, ntrees, sizeof(*locations), compare_locations_fid_file_offset);
}    


void assign_forest_ids(const int64_t ntrees, struct locations *locations, int64_t *forests, int64_t *tree_roots)
{
    /* Sort forests by tree roots -> necessary for assigning forest ids */
    sort_forests(ntrees, forests, tree_roots);
    sort_locations_on_treeroot(ntrees, locations);    
    
    /* forests and tree_roots are sorted together, on tree_roots */
    /* locations is sorted on tree roots */
    for(int64_t i=0;i<ntrees;i++) {
        XASSERT(tree_roots[i] == locations[i].tree_root,
                "tree roots[%"PRId64"] = %"PRId64" does not equal tree roots in locations = %"PRId64"\n",
                i, tree_roots[i], locations[i].tree_root);
        locations[i].forestid = forests[i];
    }
}    

struct forest_info * assign_trees_in_forest_to_same_file(const int64_t ntrees, struct locations *locations, struct locations *output_locations, const int nfiles)
{
    sort_locations_on_fid_file_offset(ntrees, locations);
    sort_locations_on_fid_file_offset(ntrees, output_locations);
    int64_t nforests = ntrees;/* Maximum value the nforests can have */

    struct forest_info *forest_info = my_calloc(sizeof(struct forest_info), 1);
    forest_info->nforests = 1;
    forest_info->forestid = my_calloc(sizeof(int64_t), nforests);
    forest_info->fileid   = my_calloc(sizeof(int64_t), nforests);
    forest_info->num_trees = my_calloc(sizeof(int64_t), nforests);
    forest_info->num_halos = my_calloc(sizeof(int64_t), nforests);
    forest_info->num_bytes = my_calloc(sizeof(int64_t), nforests);
    forest_info->offset    = my_calloc(sizeof(int64_t), nforests);
    forest_info->filename  = my_calloc(sizeof(*(forest_info->filename)), nforests);

    int64_t *histogram_fileids = my_calloc(sizeof(*histogram_fileids), nfiles);
    
    /* the fun begins here -> in case, assign all trees from a forest into the same file */
    int64_t start_forestid = locations[0].forestid;
    int64_t min_fileid = locations[0].fileid;
    int64_t max_fileid = locations[0].fileid;
    int64_t start_index_forest = 0;
    int64_t end_index_forest = 1;
    int64_t num_trees_moved = 0;

    fprintf(stderr, ANSI_COLOR_MAGENTA"Assigning all trees in a forest into the same file...."ANSI_COLOR_RESET"\n");
    /* setup the progressbar */
    int interrupted=0;
    init_my_progressbar(ntrees, &interrupted);

    int64_t forest_index = 0;
    forest_info->forestid[forest_index] = start_forestid;
    forest_info->num_trees[forest_index] = 1;
    
    for(int64_t i=1;i<ntrees;i++) {
        my_progressbar(i, &interrupted);
        if(locations[i].forestid == start_forestid) {
            nforests--;
            forest_info->num_trees[forest_index]++;
            if(locations[i].fileid < min_fileid) {
                min_fileid = locations[i].fileid;
            }
            if(locations[i].fileid > min_fileid) {
                max_fileid = locations[i].fileid;
            }
            end_index_forest++;
            continue;
        } else {
            int64_t max_common_fileid = -1;
            if(min_fileid == max_fileid) {
                for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                    output_locations[j].fileid = min_fileid;
                }
                max_common_fileid = min_fileid;
            } else {

                /* create a histogram of the fileids */
                memset(histogram_fileids, 0, sizeof(*histogram_fileids) * nfiles);
                for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                    histogram_fileids[locations[j].fileid]++;
                }

                int64_t max_common_value = 0;
                for(int j=0;j<nfiles;j++) {
                    if(histogram_fileids[j] > max_common_value) {
                        max_common_value = histogram_fileids[j];
                        max_common_fileid = j;
                    }
                }

                for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                    output_locations[j].fileid = max_common_fileid;
                    if(output_locations[j].fileid != locations[j].fileid) {
                        num_trees_moved++;
                    }
                }
            }

            /* Update the forest_info */
            XASSERT((max_common_fileid >=0 && max_common_fileid < nfiles),
                    "max common fileid = %"PRId64" is not within bounds [0,%d)\n",
                    max_common_fileid, nfiles);
                    
            forest_info->fileid[forest_index] = max_common_fileid;
            for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                forest_info->num_bytes[forest_index] += locations[j].bytes;
            }

            /* We have a new forest on the current line */
            forest_index++;
            forest_info->nforests++;
            start_forestid = locations[i].forestid;
            forest_info->num_trees[forest_index] = 1;
            forest_info->forestid[forest_index]  = start_forestid;
            start_index_forest = i;
            end_index_forest   = i+1;
            min_fileid = locations[i].fileid;
            max_fileid = locations[i].fileid;
        }
    }
    finish_myprogressbar(&interrupted);        
    free(histogram_fileids);

    if(num_trees_moved > 0) {
        fprintf(stderr,"Number of trees moved into different files = %"PRId64"\n",num_trees_moved);
    }
    fprintf(stderr, ANSI_COLOR_GREEN"Assigning all trees in a forest into the same file.......done"ANSI_COLOR_RESET"\n\n");

    XASSERT(forest_info->nforests == nforests,
            "forest_info->nforests = %"PRId64" must equal nforests = %"PRId64"\n",
            forest_info->nforests, nforests);
    
    /* free the extra memory claimed by forest_info */
    forest_info->forestid = my_realloc(forest_info->forestid, sizeof(int64_t), nforests, "struct forest_info (forestid)");
    forest_info->fileid   = my_realloc(forest_info->fileid, sizeof(int64_t), nforests, "struct forest_info (fileid)");
    forest_info->num_trees = my_realloc(forest_info->num_trees, sizeof(int64_t), nforests, "struct forest_info (num_trees)");
    forest_info->num_halos = my_realloc(forest_info->num_halos, sizeof(int64_t), nforests, "struct forest_info (num_halos)");
    forest_info->num_bytes = my_realloc(forest_info->num_bytes, sizeof(int64_t), nforests, "struct forest_info (num_bytes)");
    forest_info->offset    = my_realloc(forest_info->offset, sizeof(int64_t), nforests, "struct forest_info (offset)");
    forest_info->filename  = my_realloc(forest_info->filename, sizeof(*(forest_info->filename)), nforests, "struct forest_info (filename)");
    XASSERT(sizeof(*(forest_info->filename)) > 8,
            "size of filename holder = %zu must be larger than a pointer size",
            sizeof(*(forest_info->filename)));

    return forest_info;
}

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
                              const float inv_part_mass)
{
    int64_t nhalos = 0;
    size_t bytes_read = 0;

#ifdef USE_FGETS
    my_fseek(fp, offset, SEEK_SET);
#endif    
    
    while(bytes_read < bytes) {
        char buffer[MAXLEN+1]={'\0'};
        const int bytes_to_read = (bytes - bytes_read) > MAXLEN ? MAXLEN:(bytes - bytes_read);
#if 0        
        fprintf(stderr,"NEED TO READ %zu BYTES bytes_read = %zu bytes_to_read = %d\n", bytes, bytes_read, bytes_to_read);
#endif        

#ifdef USE_FGETS
#error FGETS does not work yet
        char *ret = fgets(buffer, bytes_to_read, fp);
        my_fread(buffer, sizeof(char), bytes_to_read, fp);
        const int64_t bytes_this_read = (ret == NULL) ? 0:(int64_t) strlen(buffer);
#else
        int64_t bytes_this_read = (int64_t) pread(fd, buffer, bytes_to_read, offset);
        if(buffer[bytes_this_read-1] != '\n') {
            for(int64_t i=bytes_this_read-2;i>=0;i--) {
                if(buffer[i] == '\n') {
                    buffer[i+1] = '\0';
                    break;
                }
                bytes_this_read = i;
            }
        }
        
        /* fprintf(stderr,"bytes_this_read = %"PRId64" len = %zu\n",bytes_this_read, strlen(buffer)); */
        /* offset += bytes_this_read; */
#endif

#if 0
        fprintf(stderr,"STARTING with NEW BUFFER. buffer = `%s'\n+++++++++++++++++++++++++++\n",buffer);
#endif        

        int64_t start_pos=0;
        while(start_pos < bytes_this_read) {
            int64_t end_pos=start_pos;
            while(end_pos < bytes_this_read && buffer[end_pos] != '\n') {
                end_pos++;
            }
            /* blank line ?*/
            if(end_pos == start_pos || (end_pos == start_pos + 1) ) {
                start_pos = end_pos + 1;
                continue;
            }

            if(end_pos < bytes_this_read) {
                XASSERT(buffer[end_pos] == '\n',"buffer[%"PRId64"] = `%c' is not a new-line\n",
                        end_pos, buffer[end_pos]);
            } 

            buffer[end_pos] = '\0';

#if 0            
            fprintf(stderr,"PARSING LINE = `%s'\n******\n",&buffer[start_pos]);
#endif
            
            /* found a new-line but before parsing check that enough memory has been allocated */
            if((forest_offset+nhalos) >= *nhalos_allocated) {
                *nhalos_allocated += 100000;
                *source_forest = my_realloc(*source_forest, sizeof(struct output_dtype), *nhalos_allocated, "struct source forest");
                *source_info = my_realloc(*source_info, sizeof(struct additional_info), *nhalos_allocated, "struct source info");
            }

            struct output_dtype *forest = *source_forest + forest_offset;
            struct additional_info *info = *source_info + forest_offset;
            const int nitems_expected = 21;
            const int nitems = sscanf(&buffer[start_pos],
                                      " %lf %"SCNd64" %lf %"SCNd64" %*d "
                                      "%"SCNd64" %"SCNd64" %*d %*d "

                                      /*the first two fields are sam_mvir and mvir. one of the two need to be assigned to Mvir.
                                        I have chosen Mvir from halofinder (second column in "%*f %f") */
                                      "%*f %f %*f %*f %f "
                                      "%*d %*f %f "
                                      "%f %f %f "
                                      "%f %f %f "
                                      "%f %f %f "
                                      "%*f %*d %*d %*d %*d %d "
                                      "%*d %*d %*d %*f %*f %f %f ",
                                      &info[nhalos].scale,
                                      &info[nhalos].id,
                                      &info[nhalos].desc_scale,
                                      &info[nhalos].descid,
                                      
                                      &info[nhalos].pid,
                                      &info[nhalos].upid,

                                      &forest[nhalos].Mvir,
                                      &forest[nhalos].VelDisp,

                                      &forest[nhalos].Vmax,

                                      &(forest[nhalos].Pos[0]), &(forest[nhalos].Pos[1]), &(forest[nhalos].Pos[2]),
                                      &(forest[nhalos].Vel[0]), &(forest[nhalos].Vel[1]), &(forest[nhalos].Vel[2]),

                                      &(forest[nhalos].Spin[0]), &(forest[nhalos].Spin[1]), &(forest[nhalos].Spin[2]),

                                      &forest[nhalos].SnapNum,

                                      &forest[nhalos].M_Mean200, &forest[nhalos].M_TopHat);
                                      
            if(nitems == nitems_expected) {
                /* correctly parsed a halo. Fix things into LHalotree convention */
                const float inv_halo_mass = 1.0f/forest[nhalos].Mvir;
                for(int k=0;k<3;k++) {
                    forest[nhalos].Spin[k] *= inv_halo_mass;
                }
                /* Convert masses to 10^10 Msun/h */
                forest[nhalos].Mvir  *= 1e-10;
                forest[nhalos].M_Mean200 *= 1e-10;
                forest[nhalos].M_TopHat *= 1e-10;

                /* Calculate the (approx.) number of particles in this halo */
                forest[nhalos].Len   = (int) roundf(forest[nhalos].Mvir * inv_part_mass);

                /* Initialize other fields to indicate they are not populated */
                forest[nhalos].FileNr = -1;
                forest[nhalos].SubhaloIndex = (int) (forest_offset + nhalos);
                forest[nhalos].SubHalfMass  = -1.0f;
                forest[nhalos].MostBoundID  = info[nhalos].id;

                /* All the mergertree indices */
                forest[nhalos].Descendant = -1;
                forest[nhalos].FirstProgenitor = -1;
                forest[nhalos].NextProgenitor = -1;
                forest[nhalos].FirstHaloInFOFgroup = -1;
                forest[nhalos].NextHaloInFOFgroup = -1;

#if 0                
                fprintf(stderr,"MY PARSING produces: "
                        "%lf %"PRId64" %lf %"PRId64" "
                        "%"PRId64" %"PRId64" "
                        "%f %f "
                        "%f "
                        "%f %f %f "
                        "%f %f %f "
                        "%f %f %f "
                        "%d "
                        "%f %f \n\n",
                        info[nhalos].scale,
                        info[nhalos].id,
                        info[nhalos].desc_scale,
                        info[nhalos].descid,
                        
                        info[nhalos].pid,
                        info[nhalos].upid,

                        forest[nhalos].Mvir,
                        forest[nhalos].VelDisp,
                        
                        forest[nhalos].Vmax,

                        forest[nhalos].Pos[0], forest[nhalos].Pos[1], forest[nhalos].Pos[2],
                        forest[nhalos].Vel[0], forest[nhalos].Vel[1], forest[nhalos].Vel[2],
                        
                        forest[nhalos].Spin[0], forest[nhalos].Spin[1], forest[nhalos].Spin[2],
                        
                        forest[nhalos].SnapNum,

                        forest[nhalos].M_Mean200, forest[nhalos].M_TopHat);
#endif 

#ifndef USE_FGETS
                offset += (end_pos - start_pos + 1);
                bytes_read += (end_pos - start_pos + 1);
#endif                 
                nhalos++;
            } else {
                fprintf(stderr,"could not parse buffer = `%s' correctly. expected = %d found = %d\n",
                        &buffer[start_pos], nitems_expected, nitems);
                fprintf(stderr,"start_pos = %"PRId64" end_pos = %"PRId64"\n",start_pos, end_pos);
                break;
            }
            start_pos = end_pos + 1;
#if 0            
            fprintf(stderr,"start_pos = %"PRId64"\n",start_pos);
#endif            
        }
    }

    return nhalos;
}    



int64_t compute_numbytes(FILE *fp, const int64_t start)
{
    const int64_t off = ftello(fp);
    return off - start; /* Or should there be a -1?*/
}    

int64_t write_forests_and_locations(const char *filename, const int64_t ntrees, const struct locations *locations)
{
    int64_t bytes = 0;
    FILE *fp = my_fopen(filename,"w");
    bytes += fprintf(fp,"############## Value-Added Consistent-Trees locations.dat file. #####################\n");
    bytes += fprintf(fp,"##    ForestID     TreeRootID     FileID    Offset        Filename          Bytes    \n");
    bytes += fprintf(fp,"#####################################################################################\n");
    
    for(int64_t i=0;i<ntrees;i++) {
        bytes += fprintf(fp," %13"PRId64" %13"PRId64" %6"PRId64" %13"PRId64"  %s %16"PRId64"\n",
                         locations[i].forestid, locations[i].tree_root, locations[i].fileid,
                         locations[i].offset, locations[i].filename, locations[i].bytes);
    }
    fclose(fp);
    return bytes;
}

int64_t find_fof_halo(const int64_t totnhalos, const struct additional_info *info, const int start_loc, const int64_t upid)
{
    int64_t loc = -1;
    while(start_loc >= 0 && start_loc < totnhalos && info[start_loc].pid != -1) {
        if(upid > info[start_loc].id) {
            for(int64_t k=start_loc+1;k<totnhalos;k++) {
                if(info[k].id == upid) {
                    loc = k;
                    break;
                }
            }
        } else {
            for(int64_t k=start_loc-1;k>=0;k--) {
                if(info[k].id == upid) {
                    loc = k;
                    break;
                }
            }
        }

        if( ! (loc >= 0 && loc < totnhalos) ) {
            return -1;
        }
        if(info[loc].pid == -1) {
            return loc;
        } else {
            return find_fof_halo(totnhalos, info, loc, info[loc].upid);
        }
    }

    return -1;
}


int fix_upid(const int64_t totnhalos, struct output_dtype *forest, struct additional_info *info, int *interrupted)
{

    int max_snapnum = -1;
    /*First sort everything on ID */
#define ID_COMPARATOR(x, y)         ((x.id > y.id ? 1:(x.id < y.id ? -1:0)))
#define SCALE_ID_COMPARATOR(x,y)    ((x.scale > y.scale ? -1:(x.scale < y.scale ? 1:ID_COMPARATOR(x, y))) )
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) {                          \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct output_dtype, forest,i,j); \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct additional_info, info, i, j) \
            }
    SGLIB_ARRAY_HEAP_SORT(struct additional_info, info, totnhalos, SCALE_ID_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);

    /* Change upid to id, so we can sort the fof's and subs to be contiguous */
    for(int64_t i=0;i<totnhalos;i++) {
        info[i].upid = (info[i].pid == -1) ? info[i].id:info[i].upid;
        if(forest[i].SnapNum > max_snapnum) {
            max_snapnum = forest[i].SnapNum;
        }
        if(info[i].pid == -1) {
            continue;
        }

        /* Only (sub)subhalos should reach here */
        /*Check if upid points to host halo with pid == -1*/
        const int64_t upid = info[i].upid;
        const int64_t loc = find_fof_halo(totnhalos, info, i, upid);
        XASSERT(loc >=0 && loc < totnhalos,
                "could not locate fof halo for i = %"PRId64" id = %"PRId64" upid = %"PRId64" loc=%"PRId64"\n",
                i, info[i].id, upid, loc);
        const int64_t new_upid = info[loc].id;
        if(new_upid != upid && CTREES_UPID_BUG == 0) {
            fprintf(stderr, ANSI_COLOR_RED "Fixing upid for i=%"PRId64" original upid =%"PRId64" new fof upid = %"PRId64 ANSI_COLOR_RESET"\n"
                    ANSI_COLOR_BLUE"Email Peter Behroozi (pbehroozi@gmail.com) about the bug in Consistent Trees code"ANSI_COLOR_RESET"\n",
                    i, upid, new_upid);
            CTREES_UPID_BUG = 1;
            *interrupted = 1;
        }
        info[i].upid = new_upid;
    }
#undef ID_COMPARATOR
#undef SCALE_ID_COMPARATOR
#undef MULTIPLE_ARRAY_EXCHANGER    

    return max_snapnum;
    
}

void assign_mergertree_indices(const int64_t totnhalos, struct output_dtype *forest, struct additional_info *info, const int max_snapnum)
{
    const int nsnapshots = max_snapnum + 1;
    float *scales = my_malloc(sizeof(*scales), nsnapshots);
    for(int i=0;i<nsnapshots;i++) {
        scales[i] = FLT_MAX;
    }
    int64_t *start_scale = my_calloc(sizeof(*start_scale), nsnapshots);
    for(int i=0;i<nsnapshots;i++) {
        start_scale[i] = -1;
    }
    int64_t *end_scale = my_calloc(sizeof(*end_scale), nsnapshots);
    
    /* Sort the trees based on scale, upid, and pid */
    /* Descending sort on scale, and then ascending sort on upid.
       The last pid sort is so that the FOF halo comes before the (sub-)subhalos. 
     */
#define PID_COMPARATOR(x, y)         ((x.pid > y.pid ? 1:(x.pid < y.pid ? -1:0)))
#define UPID_COMPARATOR(x, y)         ((x.upid > y.upid ? 1:(x.upid < y.upid ? -1:PID_COMPARATOR(x,y))))
/* Note, the negated order in the scale comparison . this ensures descending sort */    
#define SCALE_UPID_COMPARATOR(x,y)    ((x.scale > y.scale ? -1:(x.scale < y.scale ? 1:UPID_COMPARATOR(x, y))) )
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) {                          \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct output_dtype, forest,i,j); \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct additional_info, info, i, j) \
            }
    

    /* I am using heap sort rather than qsort because the forest is already somewhat sorted. For sorted data,
       qsort behaves closer to O(N^2). 
     */
    SGLIB_ARRAY_HEAP_SORT(struct additional_info, info, totnhalos, SCALE_UPID_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);
    
    /* Fix subs of subs first */
    int64_t FirstHaloInFOFgroup=-1;
    int64_t fof_id=-1;
    for(int64_t i=0;i<totnhalos;i++) {
        /* fprintf(stderr,ANSI_COLOR_RED"FirstHaloInFOFgroup = %"PRId64 ANSI_COLOR_RESET"\n",FirstHaloInFOFgroup); */
        const int snapnum = forest[i].SnapNum;
        XASSERT(snapnum >= 0 && snapnum < nsnapshots,
                "snapnum = %d is outside range [0, %d)\n",
                snapnum, nsnapshots);
        scales[snapnum] = info[i].scale;
        end_scale[snapnum] = i;
        if(start_scale[snapnum] == -1) {
            start_scale[snapnum] = i;
        }
        
        if(info[i].pid == -1) {
            forest[i].FirstHaloInFOFgroup = i;
            forest[i].NextHaloInFOFgroup = -1;
            FirstHaloInFOFgroup = i;
            fof_id = info[i].id;
            continue;
        } else {
            if(FirstHaloInFOFgroup == -1) {
                fprintf(stderr,"About to crash\n");
                for(int64_t k=0;k<totnhalos;k++) {
                    fprintf(stderr,"%03d %12.5lf %10"PRId64" %10"PRId64" %10"PRId64" %12.4g\n",
                            forest[k].SnapNum, info[k].scale, info[k].upid, info[k].pid, info[k].id, forest[k].Mvir);
                }
            }

            XASSERT(FirstHaloInFOFgroup != -1,
                    "Processing subhalos i=%"PRId64" but have not encountered FOF yet..bug\n"
                    "id = %"PRId64" pid = %"PRId64" upid = %"PRId64" snapnum = %d\n",
                    i,info[i].id, info[i].pid, info[i].upid, forest[i].SnapNum);
            
            if((info[i].pid == fof_id) || (info[i].upid == fof_id) ) {
                forest[i].FirstHaloInFOFgroup = FirstHaloInFOFgroup;
            } else {
                /* Should not reach here..I have already sorted the forest such that the FOF appears before the subs */
                fprintf(stderr,"ERROR: the sort did not place the FOF before the subs. BUG IN CTREES OR IN SORT\n");
                for(int64_t k=0;k<totnhalos;k++) {
                    fprintf(stderr,"%03d %12.5lf %10"PRId64" %10"PRId64" %10"PRId64" %12.4g\n",
                            forest[k].SnapNum, info[k].scale, info[k].upid, info[k].pid, info[k].id, forest[k].Mvir);
                }

                fprintf(stderr,"i = %"PRId64" id = %"PRId64" pid = %"PRId64" fof_id = %"PRId64" upid = %"PRId64" FirstHaloInFOFgroup = %"PRId64"\n",
                        i, info[i].id, info[i].pid, fof_id, info[i].upid, FirstHaloInFOFgroup); 
                exit(EXIT_FAILURE);
            }
            int64_t insertion_point = FirstHaloInFOFgroup;
            while(forest[insertion_point].NextHaloInFOFgroup != -1) {
                insertion_point = forest[insertion_point].NextHaloInFOFgroup;
                XASSERT(insertion_point >=0 && insertion_point < totnhalos,
                        "Inserting next halo in FOF group into invalid index. insertion_point = %"PRId64" totnhalos = %"PRId64"\n",
                        insertion_point, totnhalos);
            }
            forest[insertion_point].NextHaloInFOFgroup = i;
        }
    }

    /* for(int i=0;i<nsnapshots;i++) { */
    /*     fprintf(stderr,"scales[%04d] = %12.5lf\n", */
    /*             i,scales[i]); */
    /* } */

    
    /* Now figure out merger tree pointers. Need to set descendant, firstprogenitor and nextprogenitor.
     */
    for(int64_t i=0;i<totnhalos;i++) {
        if(info[i].descid == -1) {
            forest[i].Descendant = -1;
            continue;
        }

        int desc_snapnum = nsnapshots-1;
        const float desc_scale = info[i].desc_scale;
        const int64_t descid = info[i].descid;
        const float max_epsilon_scale = 1.0e-4;
        while((desc_snapnum >= 0) &&
              (fabsf(scales[desc_snapnum] - desc_scale) > max_epsilon_scale) ) {
            desc_snapnum--;
        }
        XASSERT(desc_snapnum >= 0 && desc_snapnum < nsnapshots && (fabsf(scales[desc_snapnum] - desc_scale) <= 1e-4),
                "Could not locate desc_snapnum. desc_snapnum = %d nsnapshots = %d \n",
                desc_snapnum, nsnapshots);

        for(int64_t k=start_scale[desc_snapnum];k<end_scale[desc_snapnum];k++) {
            if(info[k].id != descid) {
                continue;
            }
            forest[i].Descendant = k;
            
            //Now assign first progenitor + next progenitor 
            if(forest[k].FirstProgenitor == -1) {
                forest[k].FirstProgenitor = i;
                forest[i].NextProgenitor = -1;
            } else {
                /* The descendant halo already has progenitors. Figure out the correct
                   order -- should this halo be FirstProgenitor?
                 */
                const int64_t first_prog = forest[k].FirstProgenitor;
                XASSERT(first_prog >= 0 && first_prog < totnhalos,
                        "first_prog=%"PRId64" must lie within [0, %"PRId64"\n",
                        first_prog, totnhalos);
                int64_t next_prog = first_prog;
                while(forest[next_prog].NextProgenitor != -1) {
                    next_prog = forest[next_prog].NextProgenitor;
                    XASSERT(next_prog >=0 && next_prog < totnhalos,
                            "Inserting next progenitor into invalid index. insertion_point = %"PRId64" totnhalos = %"PRId64"\n",
                            next_prog, totnhalos);
                }

                /* Now figure out where to insert the current halo in the progenitor list */
                if(forest[first_prog].Mvir > forest[i].Mvir) {
                    /* The original FirstProgenitor is more massive than the current halo
                       being processed. Just append this halo index to the nextprogenitor list. 
                     */
                    forest[next_prog].NextProgenitor = i;
                    forest[i].NextProgenitor = -1;//set the next progenitor of current halo to -1. (It should already have been initialized to -1)
                } else {
                    /* Insert this halo as first progenitor and the original first progenitor next */
                    forest[k].FirstProgenitor = i;
                    forest[i].NextProgenitor = first_prog;
                }
            }
            break;
        }
    }

    free(scales);
    free(start_scale);
    free(end_scale);
    
#undef MULTIPLE_ARRAY_EXCHANGER
#undef SCALE_UPID_COMPARATOR
#undef UPID_COMPARATOR
#undef PID_COMPARATOR    

    return;
}

int64_t ** calculate_forest_info_per_file(const int64_t nfiles, int64_t *totnforests_per_file, const struct forest_info *forest_info)
{
    const int64_t nforests = forest_info->nforests;
    int64_t **nhalos_per_forest_per_file = my_malloc(sizeof(*nhalos_per_forest_per_file), nfiles);//allocate space for nfiles int64_t *'s
    for(int64_t i=0;i<nfiles;i++) {
        totnforests_per_file[i] = 0;
    }
    
    for(int64_t i=0;i<nforests;i++){
        int64_t out_fileid = forest_info->fileid[i];
        totnforests_per_file[out_fileid]++;
    }

    for(int64_t i=0;i<nfiles;i++) {
        const int64_t numforests_this_file = totnforests_per_file[i];
        nhalos_per_forest_per_file[i] = my_calloc(sizeof(**nhalos_per_forest_per_file), numforests_this_file);
    }

    return nhalos_per_forest_per_file;
}
    
        
void fix_flybys(const int64_t totnhalos, struct output_dtype *forest, struct additional_info *info)
{
#define ID_COMPARATOR(x, y)         ((x.id > y.id ? 1:(x.id < y.id ? -1:0)))
#define SCALE_ID_COMPARATOR(x,y)    ((x.scale > y.scale ? -1:(x.scale < y.scale ? 1:ID_COMPARATOR(x, y))) )
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) {                          \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct output_dtype, forest,i,j); \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct additional_info, info, i, j) \
            }
    SGLIB_ARRAY_HEAP_SORT(struct additional_info, info, totnhalos, SCALE_ID_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);

    float max_scale = info[0].scale;
    int64_t last_halo_with_max_scale = 0;
    int64_t num_fofs_last_scale = 0;
    for(int64_t i=0;i<totnhalos;i++) {
        if(info[i].scale < max_scale) {
            break;
        }
        last_halo_with_max_scale = i;
        num_fofs_last_scale += (info[i].pid == -1) ? 1:0;
    }

    XASSERT(num_fofs_last_scale > 0,
            "Weird! There must be non-zero FOF's at z=0\n. Found num_fof = %"PRId64"\n",
            num_fofs_last_scale);
    
    /* Is there anything to do? If there is only one FOF at z=0, then simply return */
    if(num_fofs_last_scale == 1) {
        return;
    }
    
    int64_t max_mass_fof_loc = 0;
    float max_mass_fof = forest[max_mass_fof_loc].Mvir;
    int64_t fof_id = info[max_mass_fof_loc].id;
    for(int64_t i=0;i<last_halo_with_max_scale;i++) {
        if(forest[i].Mvir > max_mass_fof) {
            max_mass_fof_loc = i;
            max_mass_fof = forest[max_mass_fof_loc].Mvir;
            fof_id = info[max_mass_fof_loc].id;
        }
        
    }

    XASSERT(max_mass_fof_loc < INT_MAX,
            "Most massive FOF location=%"PRId64" must be representable within INT_MAX",
            max_mass_fof_loc);

    int FirstHaloInFOFgroup = (int) max_mass_fof_loc;
    for(int64_t i=0;i<last_halo_with_max_scale;i++) {
        if(i == FirstHaloInFOFgroup) {
            continue;
        }
        if(info[i].pid == -1) {
            info[i].pid = fof_id;
#if 0
            fprintf(stderr,"Changed pid = -1 to pid = %"PRId64" for i=%"PRId64" FirstHaloInFOFgroup =%d last_halo_max_scale=%"PRId64"\n",
                    fof_id, i, FirstHaloInFOFgroup,last_halo_with_max_scale);
#endif            
        }
        info[i].upid = fof_id;
    }

#undef ID_COMPARATOR
#undef SCALE_ID_COMPARATOR
#undef MULTIPLE_ARRAY_EXCHANGER

}

int main(int argc, char **argv)
{
    char *input_dir, *output_dir;
    double part_mass = 0.0f;
    if(argc != 4) {
        usage(argc, argv);
        return EXIT_FAILURE;
    } else {
        input_dir  = argv[1];
        output_dir = argv[2];
        part_mass = atof(argv[3]);
        fprintf(stderr,"\nRunning `%s' with the following parameters \n", argv[0]);
        const int len1 = strlen(input_dir);
        const int len2 = strlen(output_dir);
        int max_str_len = len1 > len2 ? len1:len2;
        const int prefix_len = 23;
        max_str_len += prefix_len;
        const char c = '-';
        fprintf(stderr,"\t");
        for(int i=0;i<max_str_len;i++) {
            fprintf(stderr,"%c",c);
        }
        fprintf(stderr,"\n");
        fprintf(stderr,"\t Input directory   = `%s'\n",input_dir);
        fprintf(stderr,"\t Output directory  = `%s'\n",output_dir);
        fprintf(stderr,"\t Particle Mass     = `%lf'"ANSI_COLOR_BLUE" [10^10 Msun/h]"ANSI_COLOR_RESET"\n",part_mass);
        fprintf(stderr,"\t");
        for(int i=0;i<max_str_len;i++) {
            fprintf(stderr,"%c",c);
        }
        fprintf(stderr,"\n");
        fprintf(stderr,"\n");
    }
    
    if(strcmp(input_dir, output_dir) == 0) {
        fprintf(stderr,"ERROR: Input and output directories are the same..exiting\n");
        return EXIT_FAILURE;
    }
    if(part_mass > 1.0 || part_mass <= 0.0) {
        fprintf(stderr,"Expect particle mass in 10^10 Msun/h units and must be non-zero\n");
        exit(EXIT_FAILURE);
    }
    const float inv_part_mass = 1.0/part_mass;
    
    struct timeval tstart, tend;
    gettimeofday(&tstart, NULL);
    char locations_filename[MAXLEN], forests_filename[MAXLEN];
    int64_t *forests=NULL, *tree_roots=NULL;
    my_snprintf(locations_filename, MAXLEN, "%s/locations.dat", input_dir);
    my_snprintf(forests_filename, MAXLEN, "%s/forests.list", input_dir);
    fprintf(stderr, ANSI_COLOR_MAGENTA"Reading forests...."ANSI_COLOR_RESET"\n");
    const int64_t ntrees = read_forests(forests_filename, &forests, &tree_roots);
    fprintf(stderr, ANSI_COLOR_GREEN"Reading forests......done"ANSI_COLOR_RESET"\n\n");
    /* fprintf(stderr, "Number of trees = %"PRId64"\n\n",ntrees); */

    struct locations *locations = my_malloc(sizeof(*locations), ntrees);
    int nfiles = 0, BOX_DIVISIONS=0;
    fprintf(stderr, ANSI_COLOR_MAGENTA"Reading locations...."ANSI_COLOR_RESET"\n");
    const int64_t ntrees_loc = read_locations(locations_filename, ntrees, locations, &nfiles, &BOX_DIVISIONS);
    fprintf(stderr, ANSI_COLOR_GREEN"Reading locations......done"ANSI_COLOR_RESET"\n\n");
    XASSERT(ntrees == ntrees_loc,
            "ntrees=%"PRId64" should be equal to ntrees_loc=%"PRId64"\n",
            ntrees, ntrees_loc);    

    /* the following function will sort locations and forests based on tree root id*/
    assign_forest_ids(ntrees, locations, forests, tree_roots);
    
    /* Forests are now contained inside locations -> free the pointers */
    free(forests);free(tree_roots);


    FILE **tree_outputs = my_malloc(sizeof(FILE *), nfiles);
    FILE **tree_inputs  = my_malloc(sizeof(FILE *), nfiles);

    int *tree_inputs_fd = my_malloc(sizeof(*tree_inputs_fd), nfiles);
    int *tree_outputs_fd = my_malloc(sizeof(*tree_outputs_fd), nfiles);    
    char (*tree_outputs_fname)[MAXLEN]  = my_malloc(sizeof(*tree_outputs_fname), nfiles);
                                                           
    int *tree_outputs_fd_offset = my_calloc(sizeof(*tree_outputs_fd_offset), nfiles);

    int64_t *tree_counts = my_calloc(sizeof(*tree_counts), nfiles);
    int64_t *inp_file_sizes = my_calloc(sizeof(*inp_file_sizes), nfiles);
    char buffer[MAXLEN];
    for (int i=0; i<BOX_DIVISIONS; i++) {
        for (int j=0; j<BOX_DIVISIONS; j++) {
            for(int k=0; k<BOX_DIVISIONS; k++) {
                my_snprintf(buffer,MAXLEN,"%s/tree_%d_%d_%d.dat", input_dir, i, j, k);
                int id = id = i*BOX_DIVISIONS*BOX_DIVISIONS + j*BOX_DIVISIONS + k;
                tree_inputs[id]  = my_fopen(buffer, "r");
                my_fseek(tree_inputs[id],0L, SEEK_END);
                inp_file_sizes[id] = ftello(tree_inputs[id]);
                rewind(tree_inputs[id]);

                tree_inputs_fd[id]  = fileno(tree_inputs[id]);

                my_snprintf(buffer,MAXLEN,"%s/lhalotree.bin.%d", output_dir, id);
                unlink(buffer);
                my_snprintf(tree_outputs_fname[id], LOCATIONS_FILENAME_SIZE, "lhalotree.bin.%d",id);
                tree_outputs[id] = my_fopen(buffer, "w");
                tree_outputs_fd[id] = fileno(tree_outputs[id]);
                tree_outputs_fd_offset[id] = 0;
            }
        }
    }


    /* the following function will sort locations based on 1) filename 2) offsets */
    sort_locations_file_offset(ntrees, locations);

    /* holder to check later that bytes have been assigned */
    for(int64_t i=0;i<ntrees;i++) {
        locations[i].bytes = -1;/* Make sure bytes is a signed type! */
    }

    /* Create a copy of current locations */    
    struct locations *output_locations = my_malloc(sizeof(*output_locations), ntrees);
    assert(sizeof(*output_locations) == sizeof(*locations) && "locations struct is varying in size! The sky is falling!!");
    memcpy(output_locations, locations, sizeof(*locations) * ntrees);

    /* figure out the byte size for each tree */
    int64_t start = locations[0].offset;
    int64_t start_fileid = locations[0].fileid;

    /* tree_roots are 64 bit integers -> max digits in decimal = log10(2^64) < 20.
       Add 1 char for +-, in case consistent tree changes. and then strlen('#tree ')
       and the previous \n. I need to read up to previous newline.
    */
    const int64_t guess_max_linesize = 20 + 1 + 6 + 1;
    fprintf(stderr, ANSI_COLOR_MAGENTA"Calculating the number of bytes for each tree...."ANSI_COLOR_RESET"\n");
    /* setup the progressbar */
    int interrupted=0;
    init_my_progressbar(ntrees, &interrupted);

    for(int64_t i=1;i<=ntrees-1;i++) {
        my_progressbar(i, &interrupted);
        const int64_t fileid = locations[i].fileid;
        
        /* Are we starting on a new file ?*/
        if(start_fileid != fileid) {
            /* fill out the bytes for the last tree in the previous file */
            my_fseek(tree_inputs[start_fileid], 0L, SEEK_END);
            const int64_t num_bytes = compute_numbytes(tree_inputs[start_fileid], start);
            locations[i-1].bytes = num_bytes;
            output_locations[i-1].bytes = num_bytes;

            /* now we reset the start fields */
            start = locations[i].offset;
            start_fileid = locations[i].fileid;
            continue;
        }
        const int64_t current_offset_guess = locations[i].offset - guess_max_linesize;
        my_fseek(tree_inputs[fileid], current_offset_guess, SEEK_SET);
        while(1) {
            const int a = fgetc(tree_inputs[fileid]);
            if(a == EOF) {
                fprintf(stderr,"Encountered EOF while looking for end of current tree\n");
                exit(EXIT_FAILURE);
            }
            const unsigned char c = (unsigned char) a;
            if(c == '\n') {
                const int64_t num_bytes = compute_numbytes(tree_inputs[start_fileid], start);
                locations[i-1].bytes = num_bytes;
                output_locations[i-1].bytes = num_bytes;
                /* fprintf(stderr,"%"PRId64"\n",num_bytes); */
                start = locations[i].offset;
                break;
            }
        }
    }

    /* fill out the bytes for the last tree */
    {
        start = locations[ntrees-1].offset;
        const int64_t fileid = locations[ntrees-1].fileid;
        my_fseek(tree_inputs[fileid], 0L, SEEK_END);
        const int64_t num_bytes = compute_numbytes(tree_inputs[fileid], start);
        locations[ntrees-1].bytes = num_bytes;
        output_locations[ntrees-1].bytes = num_bytes;
    }
    finish_myprogressbar(&interrupted);        
    fprintf(stderr, ANSI_COLOR_GREEN"Calculating the number of bytes for each tree.....done"ANSI_COLOR_RESET"\n\n");

    /* Check that all the previous computations with locations have been copied to output_locations */
    for(int64_t i=ntrees-1;i>=0;i--) {
        XASSERT(locations[i].bytes > 0,
                "locations[%"PRId64"].bytes = %"PRId64" should be positive\n",
                i,locations[i].bytes);

        XASSERT(output_locations[i].bytes == locations[i].bytes,
                "locations[%"PRId64"].bytes = %"PRId64" should be equal output_locations->bytes = %"PRId64"\n",
                i,locations[i].bytes,output_locations[i].bytes);
        XASSERT(strncmp(output_locations[i].filename, locations[i].filename, LOCATIONS_FILENAME_SIZE) == 0,
                "output_locations[%"PRId64"].filename = %s should equal locations filename = %s\n",
                i, output_locations[i].filename, locations[i].filename);

        
        assert(output_locations[i].forestid == locations[i].forestid);
        assert(output_locations[i].tree_root == locations[i].tree_root);
        assert(output_locations[i].fileid == locations[i].fileid);
        assert(output_locations[i].offset == locations[i].offset);
        assert(output_locations[i].bytes == locations[i].bytes);
    }
    
    /* Check that the preceeding bytes computation is correct */
    {
        int64_t *total_tree_bytes = my_calloc(sizeof(*total_tree_bytes), nfiles);
        for(int64_t i=0;i<ntrees;i++) {
            /* add the number of bytes for tree in each file */
            total_tree_bytes[locations[i].fileid] += locations[i].bytes;
        }
        
        for(int i=0;i<nfiles;i++) {
            XASSERT(total_tree_bytes[i] < inp_file_sizes[i],
                    "Bytes in tree = %"PRId64" must be smaller than file size = %"PRId64"\n",
                    total_tree_bytes[i], inp_file_sizes[i]);
        }
        free(total_tree_bytes);
    }

    
    /* Now assign all trees in the same forest to the same file
       The output fileids goes into output_locations (which is otherwise a copy of locations).
     */
    struct forest_info *forest_info = assign_trees_in_forest_to_same_file(ntrees, locations, output_locations, nfiles);
    const int64_t nforests = forest_info->nforests;
    /* Fix the output filenames in output_locations */
    for(int64_t i=0;i<ntrees;i++) {
        const int64_t out_fileid = output_locations[i].fileid;
        my_snprintf(output_locations[i].filename, LOCATIONS_FILENAME_SIZE, "%s", tree_outputs_fname[out_fileid]);
    
    }
    for(int64_t i=0;i<nforests;i++) {
        const int64_t out_fileid = forest_info->fileid[i];
        my_snprintf(forest_info->filename[i], LOCATIONS_FILENAME_SIZE, "%s", tree_outputs_fname[out_fileid]);
    }
    

    int64_t *totnforests_per_file = my_calloc(sizeof(*totnforests_per_file), nfiles);
    int64_t *totnhalos_per_file   = my_calloc(sizeof(*totnhalos_per_file), nfiles);
    int64_t **nhalos_per_forest_per_file = calculate_forest_info_per_file(nfiles, totnforests_per_file, forest_info);
    int *nforests_written_per_file = my_calloc(sizeof(*nforests_written_per_file), nfiles);
    int *nhalos_written_per_file   = my_calloc(sizeof(*nhalos_written_per_file), nfiles);
    off_t *offsets_per_file        = my_calloc(sizeof(*offsets_per_file), nfiles);
    /* write the place-holders. Should I use posix_fallocate instead? */
    for(int i=0;i<nfiles;i++) {
        FILE *fp = tree_outputs[i];
        rewind(fp);
        const int zero = 0;
        my_fwrite(&zero, sizeof(int), 1, fp);//totnforests in this file 
        my_fwrite(&zero, sizeof(int), 1, fp);//totnhalos in this file 
        for(int64_t j=0;j<totnforests_per_file[i];j++) {
            my_fwrite(&zero, sizeof(int), 1, fp);//one zero for each forest in this file. will contain the number of halos in forest.
        }
        fflush(fp);
        fsync(fileno(fp));/* ensure data really are flushed to disk */
        offsets_per_file[i] = ftello(fp);
    }
    
    fprintf(stderr, ANSI_COLOR_MAGENTA"Writing out (nforests=%"PRId64") in LHALOTREE format...."ANSI_COLOR_RESET"\n", nforests);
    interrupted=0;
    init_my_progressbar(nforests, &interrupted);

    /* Now copy each one of the forests. Note a forest can have multiple trees */
    int64_t nhalos_allocated = 1000000;//1 million halos
    struct output_dtype *forest  = my_malloc(sizeof(*forest), nhalos_allocated);
    struct additional_info *info = my_malloc(sizeof(*info), nhalos_allocated);
    int64_t tree_index = 0;
    for(int64_t i=0;i<nforests;i++) {
        my_progressbar(i, &interrupted);
        int64_t forest_offset = 0;
        for(int64_t j=0;j<forest_info->num_trees[i];j++) {
            int64_t fileid = locations[tree_index].fileid;
            const int64_t nhalos = read_tree_into_forest(&nhalos_allocated, &forest, forest_offset, &info,
#ifdef USE_FGETS
                                                         tree_inputs[fileid],
#else
                                                         tree_inputs_fd[fileid],
                                                         locations[tree_index].offset,
#endif                                                         
                                                         locations[tree_index].bytes,
                                                         inv_part_mass);

            tree_index++;
            forest_offset += nhalos;
        }

        const int64_t totnhalos = forest_offset;
        const int64_t out_fileid = forest_info->fileid[i];
        const int forestindex_thisfile = nforests_written_per_file[out_fileid];
        nhalos_per_forest_per_file[out_fileid][forestindex_thisfile] = totnhalos;
        totnhalos_per_file[out_fileid] += totnhalos;
        nforests_written_per_file[out_fileid]++;
        nhalos_written_per_file[out_fileid] += totnhalos;

        /* Fix flybys -> multiple roots at z=0 must be joined such that only one root remains */
        fix_flybys(totnhalos, forest, info);

        /* Entire tree is loaded in. Fix upid's*/
        const int max_snapnum = fix_upid(totnhalos, forest, info, &interrupted);
        
        /* Now the entire tree is loaded in. Assign the mergertree indices */
        /* fprintf(stderr,"ASSIGNING MERGERTREE INDICES\n"); */

        assign_mergertree_indices(totnhalos, forest, info, max_snapnum);

        const int64_t num_bytes = sizeof(struct output_dtype) * totnhalos;
        forest_info->num_bytes[i]  = num_bytes;
        forest_info->offset[i]    = offsets_per_file[out_fileid];
        forest_info->num_halos[i]  = totnhalos;
        offsets_per_file[out_fileid] += num_bytes;
        
#if 0
        for(int64_t kk=0;kk<totnhalos;kk++) {
            XASSERT((forest[kk].Descendant == -1) || (forest[kk].Descendant >= 0 && forest[kk].Descendant < totnhalos),
                    "forest[%"PRId64"].Descendant = %d not in range\n",kk, forest[kk].Descendant);
            XASSERT((forest[kk].FirstProgenitor == -1) || (forest[kk].FirstProgenitor >= 0 && forest[kk].FirstProgenitor < totnhalos),
                    "forest[%"PRId64"].FirstProgenitor = %d not in range\n",kk, forest[kk].FirstProgenitor);
            XASSERT((forest[kk].NextProgenitor == -1) || (forest[kk].NextProgenitor >= 0 && forest[kk].NextProgenitor < totnhalos),
                    "forest[%"PRId64"].NextProgenitor = %d not in range\n",kk, forest[kk].NextProgenitor);
            XASSERT((forest[kk].FirstHaloInFOFgroup >= 0 && forest[kk].FirstHaloInFOFgroup < totnhalos),
                    "forest[%"PRId64"].FirstHaloInFOFgroup = %d not in range\n",kk, forest[kk].FirstHaloInFOFgroup);
            XASSERT((forest[kk].NextHaloInFOFgroup == -1) || (forest[kk].NextHaloInFOFgroup >= 0 && forest[kk].NextHaloInFOFgroup < totnhalos),
                    "forest[%"PRId64"].NextHaloInFOFgroup = %d not in range\n",kk, forest[kk].NextHaloInFOFgroup);

            fprintf(stderr," %12.5lf   %12d   %12"PRId64" %12"PRId64"  %12"PRId64" %14.6lf "
                    "%14.6lf %14.6lf %12.5lf %12.5lf %12.5lf "
                    "%12.4lf %12.4lf %12.4lf "
                    "%12.4lf %12.4lf %12.4lf "
                    "%12.4lf %12.4lf %12d \n",
                    info[kk].scale, forest[kk].SnapNum, info[kk].id, info[kk].pid, info[kk].upid, forest[kk].Mvir,
                    forest[kk].VelDisp, forest[kk].Vmax, forest[kk].Pos[0], forest[kk].Pos[1], forest[kk].Pos[2],
                    forest[kk].Vel[0], forest[kk].Vel[1], forest[kk].Vel[2],
                    forest[kk].Spin[0], forest[kk].Spin[1], forest[kk].Spin[2],
                    forest[kk].M_Mean200, forest[kk].M_TopHat, forest[kk].Len);
        }
#endif
        
        /* fprintf(stderr,"WRiting forest with %"PRId64" halos...\n",totnhalos); */
        my_fwrite(forest, sizeof(struct output_dtype), totnhalos, tree_outputs[out_fileid]);
        /* fprintf(stderr,"WRiting forest with %"PRId64" halos......done\n",totnhalos); */
    }
    
    for(int i=0;i<nfiles;i++) {
        FILE *fp = tree_outputs[i];
        fflush(fp);
        fsync(fileno(fp));
        rewind(fp);

        /* Check that totnforests has not overflown */
        if(totnforests_per_file[i] > INT_MAX) {
            fprintf(stderr,"\nIn file `%s' number of trees=%"PRId64" has overflown INT_MAX. Writing out garbage\n",
                    tree_outputs_fname[i], totnforests_per_file[i]);
        } else {
            if(nforests_written_per_file[i] != totnforests_per_file[i] ){
                fprintf(stderr,"\nIn file `%s' i = %d nforests_written = %d does not agree with (int64_t) totnforests_per_file = %"PRId64"\n"
                        "Something went wrong during writing the files\n",
                        tree_outputs_fname[i], i, nforests_written_per_file[i], totnforests_per_file[i]); 
            }
        }
        my_fwrite(&(nforests_written_per_file[i]), sizeof(int), 1, fp);

        if(totnhalos_per_file[i] > INT_MAX) {
            fprintf(stderr,"In file `%s' number of halos=%"PRId64" has overflown INT_MAX. Writing out garbage\n",
                    tree_outputs_fname[i], totnhalos_per_file[i]);
        } else {
            if(nhalos_written_per_file[i] != totnhalos_per_file[i]) {
                fprintf(stderr,"\nIn file `%s' i = %d nhalos_written = %d does not agree with (int64_t) totnhalos_per_file = %"PRId64"\n"
                        "Something went wrong during writing the files\n",
                        tree_outputs_fname[i], i, nhalos_written_per_file[i], totnhalos_per_file[i]); 
            }
        }

        my_fwrite(&(nhalos_written_per_file[i]), sizeof(int), 1, fp);
        const int64_t *nhalos_this_file = nhalos_per_forest_per_file[i];
        for(int64_t j=0;j<totnforests_per_file[i];j++) {
            if(nhalos_this_file[j] > INT_MAX) {
                fprintf(stderr,"\nIn file `%s' number of halos=%"PRId64" in tree=%"PRId64" has overflown INT_MAX. Writing out garbage\n",
                        tree_outputs_fname[i], nhalos_this_file[j], j);
            }
            const int nhalos_this_tree = (int) nhalos_this_file[j];
            my_fwrite(&nhalos_this_tree, sizeof(nhalos_this_tree), 1, fp);
        }
    }
    finish_myprogressbar(&interrupted);
    fprintf(stderr, ANSI_COLOR_GREEN"Writing out (nforests=%"PRId64") in LHALOTREE format.......done"ANSI_COLOR_RESET"\n", nforests);

    /* Write out a new offsets file that can serve to connect back the LHaloTree files to the original consistent tree files */
    {
        my_snprintf(buffer, MAXLEN, "%s/lhalotree_offsets.dat", output_dir);
        FILE *fp = my_fopen(buffer,"w");
        fprintf(fp,"# ForestID   FileID   Filename   Offset  NumBytes   NumHalos \n");
        fprintf(fp,"##############################################################################################################\n");
        fprintf(fp,"#         ForestID   FileID       Filename                  Offset               NumBytes          NumHalos   \n");
        fprintf(fp,"##############################################################################################################\n");

        for(int64_t i=0;i<forest_info->nforests;i++) {
            fprintf(fp,"%18"PRId64" %10"PRId64"  %s   %18"PRId64"  %18"PRId64"  %18"PRId64"\n",
                    forest_info->forestid[i],forest_info->fileid[i], forest_info->filename[i],
                    forest_info->offset[i], forest_info->num_bytes[i], forest_info->num_halos[i]);
        }
        fclose(fp);
    }

    

    /* close open file pointers + free memory for file pointers */
    for(int i=0;i<nfiles;i++) {
        fclose(tree_inputs[i]);
        fclose(tree_outputs[i]);
    }
    free(tree_inputs);free(tree_outputs);
    free(tree_inputs_fd);free(tree_outputs_fd);

    /* free other heap allocations */
    free(tree_outputs_fd_offset);
    free(tree_outputs_fname);
    free(tree_counts);
    free(inp_file_sizes);
    free(locations);
    free(output_locations);

    free(forest_info->forestid);
    free(forest_info->fileid);
    free(forest_info->num_trees);
    free(forest_info->num_halos);
    free(forest_info->num_bytes);
    free(forest_info->offset);
    free(forest_info->filename);
    free(forest_info);

    for(int64_t i=0;i<nfiles;i++) {
        free(nhalos_per_forest_per_file[i]);
    }
    free(nhalos_per_forest_per_file);
    free(nforests_written_per_file);
    free(nhalos_written_per_file);
    free(offsets_per_file);
    free(totnforests_per_file);
    free(totnhalos_per_file);
    free(forest);free(info);
    
    gettimeofday(&tend, NULL);
    fprintf(stderr,"\n\nWrote out %"PRId64" trees in %"PRId64" forests in the LHALOTree format. Time taken = %0.2g seconds\n\n",
            ntrees, nforests, ADD_DIFF_TIME(tstart, tend));
    
    return EXIT_SUCCESS;
}
    
    
