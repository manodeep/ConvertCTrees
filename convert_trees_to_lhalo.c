#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "convert_trees_to_lhalo.h"
#include "sglib.h"
#include "progressbar.h"

#ifdef USE_STRINGPARSE
#include "stringparse.h"
#include "check_syscalls.h"
#endif


/* #define FIELD_CHECKER(field, Nmax, fieldname)   (XASSERT((field == -1) || (field >=0 && field < Nmax), fieldname " = %"PRId64" must be " */

int CTREES_UPID_FEATURE = 0;

int64_t read_forests(const char *filename, const char *output_dir, int64_t **f, int64_t **t)
{
    char forest_bin_file[MAXLEN];
    my_snprintf(forest_bin_file, MAXLEN, "%s/forests.bin", output_dir);
    FILE *forests_fp = fopen(forest_bin_file,"r");
    int64_t ntrees;
    if(forests_fp == NULL) {
        char buffer[MAXBUFSIZE];
        const char comment = '#';
        /* By passing the comment character, getnumlines
           will return the actual number of lines, ignoring
           the first header line. 
        */
        
        ntrees = getnumlines(filename, comment);
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
        /* Output the forests in binary */
        forests_fp = my_fopen(forest_bin_file,"w");
        size_t dummy = sizeof(*forests);
        my_fwrite(&dummy, sizeof(dummy), 1, forests_fp);
        my_fwrite(&ntrees, sizeof(ntrees), 1, forests_fp);
        my_fwrite(forests, sizeof(*forests), ntrees, forests_fp);
        my_fwrite(tree_roots, sizeof(*tree_roots), ntrees, forests_fp);
        fclose(forests_fp);
    } else {
        /*Found the binary forests file -> read it in */
        size_t dummy;
        my_fread(&dummy, sizeof(dummy), 1, forests_fp);
        my_fread(&ntrees, sizeof(ntrees), 1, forests_fp);
        assert(ntrees >= 0);
        *f = my_malloc(sizeof(int64_t), ntrees);
        *t = my_malloc(sizeof(int64_t), ntrees);
        int64_t *forests    = *f;
        int64_t *tree_roots = *t;
        assert(dummy == sizeof(*forests));

        my_fread(forests, sizeof(*forests), ntrees, forests_fp);
        my_fread(tree_roots, sizeof(*tree_roots), ntrees, forests_fp);
        fclose(forests_fp);
    }
        
    return ntrees;
}    


int64_t read_locations(const char *filename, const int64_t ntrees, struct locations *l, int64_t *num_files, int64_t *box_div)
{
    char buffer[MAXBUFSIZE];
    int64_t max_fileid = 0;
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

    for(int64_t i=0;i<ntrees_found;i++){
        if (locations[i].fileid > max_fileid) {
            max_fileid = locations[i].fileid;
        }
    }

    /* number of files is one greater from 0 based indexing of C files */
    *num_files = max_fileid + 1;
    const int box_divisions = (int) round(cbrt(*num_files));
    const int box_cube = box_divisions * box_divisions * box_divisions;
    XASSERT( (box_cube) == (*num_files),
             "box_divisions^3=%d should be equal to nfiles=%"PRId64"\n",
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

struct forest_info * assign_trees_in_forest_to_same_file(const int64_t ntrees, struct locations *locations, struct locations *output_locations, const int64_t nfiles)
{
    sort_locations_on_fid_file_offset(ntrees, locations);
    sort_locations_on_fid_file_offset(ntrees, output_locations);
    int64_t nforests = ntrees;/* Maximum value the nforests can have */

    struct forest_info *forest_info = my_calloc(sizeof(struct forest_info), 1);
    forest_info->nforests = 1;
    forest_info->forestid = my_calloc(sizeof(*(forest_info->forestid)), nforests);
    forest_info->fileid   = my_calloc(sizeof(*(forest_info->fileid)), nforests);
    forest_info->num_trees = my_calloc(sizeof(*(forest_info->num_trees)), nforests);
    forest_info->num_halos = my_calloc(sizeof(*(forest_info->num_halos)), nforests);
    forest_info->num_ascii_bytes = my_calloc(sizeof(*(forest_info->num_ascii_bytes)), nforests);
    forest_info->num_binary_bytes = my_calloc(sizeof(*(forest_info->num_binary_bytes)), nforests);
    forest_info->offset    = my_calloc(sizeof(*(forest_info->offset)), nforests);
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
                /* fprintf(stderr,"For forest id = %"PRId64" trees are stored in separate files (min, max) = (%"PRId64", %"PRId64")\n", */
                /*         start_forestid, min_fileid, max_fileid); */
                /* interrupted=1; */
                
                /* create a histogram of the fileids */
                memset(histogram_fileids, 0, sizeof(*histogram_fileids) * nfiles);
                for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                    histogram_fileids[locations[j].fileid]++;
                }

                int64_t max_common_value = 0;
                for(int64_t j=0;j<nfiles;j++) {
                    if(histogram_fileids[j] > max_common_value) {
                        max_common_value = histogram_fileids[j];
                        max_common_fileid = j;
                    }
                }

                for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                    output_locations[j].fileid = max_common_fileid;
                    if(output_locations[j].fileid != locations[j].fileid) {
                        num_trees_moved++;
                        /* fprintf(stderr,"Moved tree = %10"PRId64" from fileid=%3"PRId64" to fileid=%3"PRId64"\n",locations[j].tree_root, locations[j].fileid, output_locations[j].fileid); */
                        /* interrupted=1; */
                    }
                }
            }

            /* Update the forest_info */
            XASSERT((max_common_fileid >=0 && max_common_fileid < nfiles),
                    "max common fileid = %"PRId64" is not within bounds [0,%"PRId64")\n",
                    max_common_fileid, nfiles);
                    
            forest_info->fileid[forest_index] = max_common_fileid;
            for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                forest_info->num_ascii_bytes[forest_index] += locations[j].bytes;
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
    forest_info->forestid = my_realloc(forest_info->forestid, sizeof(*(forest_info->forestid)), nforests, "struct forest_info (forestid)");
    forest_info->fileid   = my_realloc(forest_info->fileid, sizeof(*(forest_info->fileid)), nforests, "struct forest_info (fileid)");
    forest_info->num_trees = my_realloc(forest_info->num_trees, sizeof(*(forest_info->num_trees)), nforests, "struct forest_info (num_trees)");
    forest_info->num_halos = my_realloc(forest_info->num_halos, sizeof(*(forest_info->num_halos)), nforests, "struct forest_info (num_halos)");
    forest_info->num_ascii_bytes = my_realloc(forest_info->num_ascii_bytes, sizeof(*(forest_info->num_ascii_bytes)), nforests, "struct forest_info (num_ascii_bytes)");
    forest_info->num_binary_bytes = my_realloc(forest_info->num_binary_bytes, sizeof(*(forest_info->num_binary_bytes)), nforests, "struct forest_info (num_binary_bytes)");
    forest_info->offset    = my_realloc(forest_info->offset, sizeof(*(forest_info->offset)), nforests, "struct forest_info (offset)");
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
        const int64_t bytes_to_read = (bytes - bytes_read) > MAXLEN ? MAXLEN:(bytes - bytes_read);
#if 0        
        fprintf(stderr,"NEED TO READ %zu BYTES bytes_read = %zu bytes_to_read = %"PRId64"\n", bytes, bytes_read, bytes_to_read);
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
                *source_forest = my_realloc(*source_forest, sizeof(**source_forest), *nhalos_allocated, "struct source forest");
                *source_info = my_realloc(*source_info, sizeof(**source_info), *nhalos_allocated, "struct source info");
            }

            struct output_dtype *forest = *source_forest + forest_offset;
            struct additional_info *info = *source_info + forest_offset;
            const int nitems_expected = 21;

#ifndef USE_STRINGPARSE
            const int nitems = sscanf(&buffer[start_pos],
                                      "%lf %"SCNd64" %lf %"SCNd64" %*d "
                                      "%"SCNd64" %"SCNd64" %*d %*d "
                                      
                                      /*the first two fields are sam_mvir and mvir. one of the two need to be assigned to Mvir.
                                        I have chosen Mvir from halofinder (second column in "%*f %f") */
                                      "%*f %f %*f %*f %f "
                                          "%*d %*f %f "
                                          "%f %f %f "
                                      "%f %f %f "
                                      "%f %f %f "
                                      "%*f %*d %*d %*d %*d %d "
                                      "%*d %*d %*f %*f %*f %f %f ",
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
#else
            buffer[end_pos] = '\n';
            //Use Peter Behroozi's custom string parser.
            SHORT_PARSETYPE;

#define NUM_INPUTS 39
            enum short_parsetype stypes[NUM_INPUTS] =  {
                LF, LD, LF, LD, K,
                LD, LD, K, K,
                K, F, K, K, F,
                K, K, F,
                F, F, F,
                F, F, F,
                F, F, F,
                K, K, K, K, K, D,
                K, K, K, K, K, F, F,
            };
            enum parsetype types[NUM_INPUTS];
            for (int ii=0; ii<NUM_INPUTS; ii++) types[ii] = stypes[ii];
            
            void *data[NUM_INPUTS] = {&(info[nhalos].scale),
                                      &(info[nhalos].id),
                                      &(info[nhalos].desc_scale),
                                      &(info[nhalos].descid),
                                      NULL,
                                      
                                      &(info[nhalos].pid),
                                      &(info[nhalos].upid),
                                      NULL,
                                      NULL,
                                      
                                      NULL,
                                      &(forest[nhalos].Mvir),
                                      NULL,
                                      NULL,
                                      &(forest[nhalos].VelDisp),
                                      
                                      NULL,
                                      NULL,
                                      &(forest[nhalos].Vmax),
                                      
                                      &(forest[nhalos].Pos[0]), &(forest[nhalos].Pos[1]), &(forest[nhalos].Pos[2]),
                                      
                                      &(forest[nhalos].Vel[0]), &(forest[nhalos].Vel[1]), &(forest[nhalos].Vel[2]),
                                      
                                      &(forest[nhalos].Spin[0]), &(forest[nhalos].Spin[1]), &(forest[nhalos].Spin[2]),
                                      
                                      NULL,
                                      NULL,
                                      NULL,
                                      NULL,
                                      NULL,
                                      &(forest[nhalos].SnapNum),
                                      
                                      NULL,
                                      NULL,
                                      NULL,
                                      NULL,
                                      NULL, 
                                      &(forest[nhalos].M_Mean200),
                                      &(forest[nhalos].M_TopHat),
            };
            //parse string

            //stringparse returns int64_t. to be consistent with the previous
            //section I made nitems as an int
            int nitems = (int) stringparse(&buffer[start_pos], data, (enum parsetype *)types, NUM_INPUTS);
            if (nitems == NUM_INPUTS) {
                nitems = nitems_expected;
            }
            buffer[end_pos] = '\0';
#undef NUM_INPUTS
            
            

#endif            

            XASSERT(nitems == nitems_expected,
                    ANSI_COLOR_RED"could not parse buffer = `%s' correctly. expected = %d found = %d\n"
                    "start_pos = %"PRId64" end_pos = %"PRId64"\n"
                    "Try deleting the locations.bin and forests.bin files in the output directory and re-running the code\n",
                    &buffer[start_pos], nitems_expected, nitems,start_pos, end_pos);
            
                                      
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

            /* Carry the Rockstar/Ctrees generated haloID through */
            forest[nhalos].MostBoundID  = info[nhalos].id;
            
            /* All the mergertree indices */
            forest[nhalos].Descendant = -1;
            forest[nhalos].FirstProgenitor = -1;
            forest[nhalos].NextProgenitor = -1;
            forest[nhalos].FirstHaloInFOFgroup = -1;
            forest[nhalos].NextHaloInFOFgroup = -1;
            /* if(info[nhalos].id == 3058982278 || info[nhalos].id == 3058982280) { */
            /*     fprintf(stderr,"info[%"PRId64"].id = %"PRId64". pid = %"PRId64" upid = %"PRId64"\n", */
            /*             nhalos, info[nhalos].id, info[nhalos].pid, info[nhalos].upid); */
            /*     fprintf(stderr,"##buffer = `%s'###\n",&buffer[start_pos]); */
            /*     fprintf(stderr,"M200 = %lf Mtophat = %lf\n", forest[nhalos].M_Mean200, forest[nhalos].M_TopHat); */
            /* } */
            

#if 0
            if(num_lines_printed < 10000) {
                fprintf(stdout,"MY PARSING produces: "
                        "%lf %"PRId64" %lf %"PRId64" "
                        "%"PRId64" %"PRId64" "
                        "%f %f "
                        "%f "
                        "%f %f %f "
                        "%f %f %f "
                        "%f %f %f "
                        "%d "
                        "%f %f \n",
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
                num_lines_printed++;
            }
#endif 
                    
#ifndef USE_FGETS
                offset += (end_pos - start_pos + 1);
                bytes_read += (end_pos - start_pos + 1);
#endif                 
                nhalos++;
            start_pos = end_pos + 1;
#if 0            
            fprintf(stderr,"start_pos = %"PRId64"\n",start_pos);
#endif            
        }
    }

    return nhalos;
}    


int64_t compute_numbytes_with_off(const int64_t off, const int64_t start)
{
    return off - start; /* Or should there be a -1?*/
}    


int64_t compute_numbytes(FILE *fp, const int64_t start)
{
    const int64_t off = ftello(fp);
    return compute_numbytes_with_off(off, start);
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

int64_t find_fof_halo(const int64_t totnhalos, const struct additional_info *info, const int start_loc, const int64_t upid, int verbose, int64_t calldepth)
{
    XASSERT(totnhalos < INT_MAX,
            "Totnhalos must be less than %d. Otherwise indexing with int (start_loc) will break\n", INT_MAX);
    int64_t loc = -1;
    assert(info[start_loc].pid != -1);
    if(calldepth >= 3) {
        verbose = 1;
    }
    assert(calldepth <= 30);
    while(start_loc >= 0 && start_loc < totnhalos && info[start_loc].pid != -1) {

        /* if(info[start_loc].id == 3058456321) { */
        /*     fprintf(stderr,"info[%d].id = 3058456321. pid = %"PRId64" upid = %"PRId64"\n", */
        /*             start_loc, info[start_loc].pid, info[start_loc].upid); */
        /* } */


        /* if(verbose == 1) { */
        /*     fprintf(stderr,"start_loc = %d id = %"PRId64" pid = %"PRId64"\n", start_loc, info[start_loc].id, info[start_loc].pid); */
        /*     fprintf(stderr,"scale = %lf pid = %"PRId64" upid = %"PRId64"\n", */
        /*             info[start_loc].scale, info[start_loc].pid, info[start_loc].upid); */
        /* } */
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
            /* if(verbose == 1) { */
            /*     fprintf(stderr,"calling find_fof_halo again with loc =%"PRId64" (int) loc = %d start_loc was =%d\n",loc, (int) loc, start_loc); */
            /*     fprintf(stderr,"scale = %lf id = %"PRId64" pid = %"PRId64" upid = %"PRId64" calldepth=%"PRId64"\n", */
            /*             info[loc].scale, info[loc].id, info[loc].pid, info[loc].upid,calldepth); */
            /* } */
            calldepth++;
            return find_fof_halo(totnhalos, info, (int) loc, info[loc].upid, verbose, calldepth);
        }
    }

    return -1;
}


int fix_upid(const int64_t totnhalos, struct output_dtype *forest, struct additional_info *info, int *interrupted, const int verbose)
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
    //I am paranoid -> so I am going to set all FOF upid's first and then
    //use the upid again. Two loops are required but that relaxes any assumptions
    //about ordering of fof/subhalos. 
    for(int64_t i=0;i<totnhalos;i++) {
        info[i].upid = (info[i].pid == -1) ? info[i].id:info[i].upid;
        if(forest[i].SnapNum > max_snapnum) {
            max_snapnum = forest[i].SnapNum;
        }
    }

    for(int64_t i=0;i<totnhalos;i++) {
        if(info[i].pid == -1) {
            continue;
        }
        
        /* Only (sub)subhalos should reach here */
        /*Check if upid points to host halo with pid == -1*/
        const int64_t upid = info[i].upid;
        int64_t calldepth=0;
        /* fprintf(stderr,"CALLING FIND FOF HALO with i = %"PRId64" id = %"PRId64" upid = %"PRId64"\n", i, info[i].id, upid); */
        const int64_t loc = find_fof_halo(totnhalos, info, i, upid, verbose, calldepth);
        XASSERT(loc >=0 && loc < totnhalos,
                "could not locate fof halo for i = %"PRId64" id = %"PRId64" upid = %"PRId64" loc=%"PRId64"\n",
                i, info[i].id, upid, loc);
        const int64_t new_upid = info[loc].id;
        if(new_upid != upid && CTREES_UPID_FEATURE == 0) {
            fprintf(stderr, ANSI_COLOR_RED "Fixing upid for i=%"PRId64" original upid =%"PRId64" new fof upid = %"PRId64 ANSI_COLOR_RESET"\n",
                    i, upid, new_upid);
            CTREES_UPID_FEATURE = 1;
            *interrupted = 1;
        }
        info[i].upid = new_upid;
        info[i].pid  = new_upid;
    }
#undef ID_COMPARATOR
#undef SCALE_ID_COMPARATOR
#undef MULTIPLE_ARRAY_EXCHANGER    

    return max_snapnum;
    
}

void assign_mergertree_indices(const int64_t totnhalos, struct output_dtype *forest, struct additional_info *info, const int max_snapnum)
{
    /* fprintf(stderr,"IN MERGERTREE totnhalos = %"PRId64"\n",totnhalos); */
    
    const int nsnapshots = max_snapnum + 1;
    double *scales = my_malloc(sizeof(*scales), nsnapshots);
    for(int i=0;i<nsnapshots;i++) {
        scales[i] = DBL_MAX;
    }
    int64_t *start_scale = my_calloc(sizeof(*start_scale), nsnapshots);
    for(int i=0;i<nsnapshots;i++) {
        start_scale[i] = -1;
    }
    int64_t *end_scale = my_calloc(sizeof(*end_scale), nsnapshots);
    
    /* Sort the trees based on scale, upid, and pid */
    /* Descending sort on scale, and then ascending sort on upid.
       The pid sort is so that the FOF halo comes before the (sub-)subhalos.
       The last id sort is such that the ordering of (sub-)subhalos
       is unique (stable sort, since id's are unique) 
     */
#define ID_COMPARATOR(x, y)         ((x.id > y.id ? 1:(x.id < y.id ? -1: 0)))    
#define PID_COMPARATOR(x, y)         ((x.pid > y.pid ? 1:(x.pid < y.pid ? -1:ID_COMPARATOR(x,y))))
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
            XASSERT(i < INT_MAX,
                    "Assigning to integer i = %"PRId64" is more than %d\n",
                    i, INT_MAX);
            forest[i].FirstHaloInFOFgroup = (int) i;
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
            
            if(info[i].upid == fof_id) {
                XASSERT(FirstHaloInFOFgroup < INT_MAX,
                        "Assigning FirstHaloInFOFgroup = %"PRId64". Must be less than %d\n",
                        FirstHaloInFOFgroup, INT_MAX);
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
                const int64_t nexthalo = forest[insertion_point].NextHaloInFOFgroup;
                XASSERT(nexthalo >=0 && nexthalo < totnhalos,
                        "Inserting next halo in FOF group into invalid index. nexthalo = %"PRId64" totnhalos = %"PRId64"\n",
                        nexthalo, totnhalos);

                /* if(forest[nexthalo].Mvir < forest[i].Mvir) { */
                /*     forest[i].NextHaloInFOFgroup = nexthalo; */
                /*     break; */
                /* } */
                insertion_point = nexthalo;
            }
            XASSERT(i < INT_MAX,
                    "Assigning FirstHaloInFOFgroup = %"PRId64". Must be less than %d\n",
                    i, INT_MAX);

            forest[insertion_point].NextHaloInFOFgroup = i;
        }
    }

    /* Now figure out merger tree pointers. Need to set descendant, firstprogenitor and nextprogenitor.
     */
    for(int64_t i=0;i<totnhalos;i++) {
        if(info[i].descid == -1) {
            forest[i].Descendant = -1;
            continue;
        }

        int desc_snapnum = nsnapshots-1;
        const double desc_scale = info[i].desc_scale;
        const int64_t descid = info[i].descid;
        const double max_epsilon_scale = 1.0e-4;
        while((desc_snapnum >= 0) &&
              (fabs(scales[desc_snapnum] - desc_scale) > max_epsilon_scale) ) {
            desc_snapnum--;
        }
        XASSERT(desc_snapnum >= 0 && desc_snapnum < nsnapshots && (fabs(scales[desc_snapnum] - desc_scale) <= 1e-4),
                "Could not locate desc_snapnum. desc_snapnum = %d nsnapshots = %d \n",
                desc_snapnum, nsnapshots);

        /*start_scale and end_scale are inclusive. Hence the stopping condition is "<=" rather than simply "<" */
        int64_t desc_loc = start_scale[desc_snapnum];
        while(desc_loc >= start_scale[desc_snapnum] && desc_loc <= end_scale[desc_snapnum] && info[desc_loc].id != descid) {
            desc_loc++;
        }
        XASSERT(desc_loc >= start_scale[desc_snapnum] && desc_loc <= end_scale[desc_snapnum],
                "Desc loc = %"PRId64" for snapnum = %d is outside range [%"PRId64", %"PRId64"]\n",
                desc_loc, desc_snapnum, start_scale[desc_snapnum], end_scale[desc_snapnum]);
        XASSERT(info[desc_loc].id == descid,
                "Should have found descendant id = %"PRId64" but info[%"PRId64"]=%"PRId64" instead \n",
                descid, desc_loc, info[desc_loc].id);

        XASSERT(desc_loc < INT_MAX,
                "desc_loc = %"PRId64" must be less than INT_MAX = %d\n",
                desc_loc, INT_MAX);
        
        forest[i].Descendant = desc_loc;
        
        //Now assign first progenitor + next progenitor 
        if(forest[desc_loc].FirstProgenitor == -1) {
            forest[desc_loc].FirstProgenitor = i;
            forest[i].NextProgenitor = -1;

        } else {
            /* The descendant halo already has progenitors. Figure out the correct
               order -- should this halo be FirstProgenitor?
               Not necessary but ensure nextprog are ordered by mass.
            */
            const int first_prog = forest[desc_loc].FirstProgenitor;
            XASSERT(first_prog >= 0 && first_prog < totnhalos,
                    "first_prog=%d must lie within [0, %"PRId64"\n",
                    first_prog, totnhalos);
            if(forest[first_prog].Mvir < forest[i].Mvir) {
                XASSERT(i < INT_MAX,
                        "Assigning Nextprogenitor = %"PRId64" to an int will result in garbage. INT_MAX = %d\n",
                        i, INT_MAX);
                forest[desc_loc].FirstProgenitor = i;
                forest[i].NextProgenitor = first_prog;
            } else {
                int64_t insertion_point = first_prog;
                while(forest[insertion_point].NextProgenitor != -1) {
                    const int64_t next_prog = forest[insertion_point].NextProgenitor;
                    XASSERT(next_prog >=0 && next_prog < totnhalos,
                            "Inserting next progenitor into invalid index. insertion_point = %"PRId64" totnhalos = %"PRId64"\n",
                            next_prog, totnhalos);
                    
                    /* if(forest[next_prog].Mvir < forest[i].Mvir) { */
                    /*     forest[i].NextProgenitor = next_prog; */
                    /*     break; */
                    /* } */
                    insertion_point = next_prog;
                }
                XASSERT(i < INT_MAX,
                        "Assigning Nextprogenitor = %"PRId64" to an int will result in garbage. INT_MAX = %d\n",
                        i, INT_MAX);
                forest[insertion_point].NextProgenitor = i;
            }
        }
    }

    free(scales);
    free(start_scale);
    free(end_scale);
    
#undef MULTIPLE_ARRAY_EXCHANGER
#undef SCALE_UPID_COMPARATOR
#undef UPID_COMPARATOR
#undef PID_COMPARATOR    
#undef ID_COMPARATOR
    
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
    
        
int fix_flybys(const int64_t totnhalos, struct output_dtype *forest, struct additional_info *info, int verbose)
{
#define ID_COMPARATOR(x, y)         ((x.id > y.id ? 1:(x.id < y.id ? -1:0)))
#define SCALE_ID_COMPARATOR(x,y)    ((x.scale > y.scale ? -1:(x.scale < y.scale ? 1:ID_COMPARATOR(x, y))) )
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) {                          \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct output_dtype, forest,i,j); \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct additional_info, info, i, j) \
            }
    SGLIB_ARRAY_HEAP_SORT(struct additional_info, info, totnhalos, SCALE_ID_COMPARATOR, MULTIPLE_ARRAY_EXCHANGER);

#undef ID_COMPARATOR
#undef SCALE_ID_COMPARATOR
#undef MULTIPLE_ARRAY_EXCHANGER    

    double max_scale = info[0].scale;
    int64_t last_halo_with_max_scale = 1;
    int64_t num_fofs_last_scale = info[0].pid == -1 ? 1:0;
    for(int64_t i=1;i<totnhalos;i++) {
        if(info[i].scale < max_scale) {
            break;
        }
        num_fofs_last_scale += (info[i].pid == -1) ? 1:0;
        last_halo_with_max_scale = i;
    }
    if(num_fofs_last_scale == 0) {
        fprintf(stderr,"ERROR: NO FOFs at max scale = %lf Will crash - here's some info that might help debug\n", max_scale);
        fprintf(stderr,"Last scale halo id (likely tree root id ) = %"PRId64" at a = %lf\n",info[0].id, info[0].scale);
        fprintf(stderr,"########################################################\n");
        fprintf(stderr,"# snap     id      pid      upid    mass     scale      \n");
        fprintf(stderr,"########################################################\n");
        for(int64_t i=0;i<=last_halo_with_max_scale;i++) {
            fprintf(stderr,"%d  %10"PRId64"  %10"PRId64" %10"PRId64" %12.6e  %20.8e\n",
                    forest[i].SnapNum, info[i].id, info[i].pid, info[i].upid, forest[i].Mvir, info[i].scale);
        }
        fprintf(stderr,"All halos now:\n\n");
        for(int64_t i=0;i<totnhalos;i++) {
            fprintf(stderr,"%d  %10"PRId64"  %10"PRId64" %10"PRId64" %12.6e %20.8e\n",
                    forest[i].SnapNum, info[i].id, info[i].pid, info[i].upid, forest[i].Mvir, info[i].scale);
        }
        return -1;
    }
    
    /* Is there anything to do? If there is only one FOF at z=0, then simply return */
    if(num_fofs_last_scale == 1) {
        return EXIT_SUCCESS;
    }
    
    int64_t max_mass_fof_loc = -1;
    float max_mass_fof = -1.0f;
    int64_t fof_id = -1;
    for(int64_t i=0;i<=last_halo_with_max_scale;i++) {
        if(forest[i].Mvir > max_mass_fof && info[i].pid == -1) {
            max_mass_fof_loc = i;
            max_mass_fof = forest[max_mass_fof_loc].Mvir;
            fof_id = info[max_mass_fof_loc].id;
        }
    }

    XASSERT(fof_id != -1,
            "There must be at least one FOF halo.");
    XASSERT(max_mass_fof_loc < INT_MAX,
            "Most massive FOF location=%"PRId64" must be representable within INT_MAX=%d",
            max_mass_fof_loc, INT_MAX);

    int FirstHaloInFOFgroup = (int) max_mass_fof_loc;
    for(int64_t i=0;i<=last_halo_with_max_scale;i++) {
        if(i == FirstHaloInFOFgroup) {
            continue;
        }
        if(info[i].pid == -1) {
            //Show that this halo was switched from being a central
            //just flip the sign. (MostBoundID should not have negative
            //values -> this would signify a flyby)
            forest[i].MostBoundID = -forest[i].MostBoundID;
            info[i].pid = fof_id;
            if(verbose == 1) {
                fprintf(stderr,"id = %"PRId64" changed pid = -1 to pid = %"PRId64" for i=%"PRId64" FirstHaloInFOFgroup =%d last_halo_max_scale=%"PRId64"\n",
                        info[i].id, fof_id, i, FirstHaloInFOFgroup,last_halo_with_max_scale);
            }
        }
        info[i].upid = fof_id;
    }

    return EXIT_SUCCESS;
}

void validate_fields(const int64_t totnhalos, const struct output_dtype *forest, const struct additional_info *info, const int max_snapnum)
{

    FILE *fp = my_fopen("mergertree","a+");
    fprintf(fp,"## Index  Snap   Desc    FiProg  neprog  Fof   nexthalo  Mvir     id\n");

    for(int64_t i=0;i<totnhalos;i++) {
        XASSERT(((forest[i].Descendant == -1 && info[i].descid == -1) || (forest[i].Descendant >= 0 && forest[i].Descendant < totnhalos)),
                "forest[%"PRId64"].Descendant = %d not in range\n",i, forest[i].Descendant);
        XASSERT(((forest[i].FirstProgenitor == -1) || (forest[i].FirstProgenitor >= 0 && forest[i].FirstProgenitor < totnhalos)),
                "forest[%"PRId64"].FirstProgenitor = %d not in range\n",i, forest[i].FirstProgenitor);
        XASSERT(((forest[i].NextProgenitor == -1) || (forest[i].NextProgenitor >= 0 && forest[i].NextProgenitor < totnhalos)),
                "forest[%"PRId64"].NextProgenitor = %d not in range\n",i, forest[i].NextProgenitor);
        XASSERT(((forest[i].FirstHaloInFOFgroup >= 0 && forest[i].FirstHaloInFOFgroup < totnhalos)),
                "forest[%"PRId64"].FirstHaloInFOFgroup = %d not in range\n",i, forest[i].FirstHaloInFOFgroup);
        XASSERT(((forest[i].NextHaloInFOFgroup == -1) || (forest[i].NextHaloInFOFgroup >= 0 && forest[i].NextHaloInFOFgroup < totnhalos)),
                "forest[%"PRId64"].NextHaloInFOFgroup = %d not in range\n",i, forest[i].NextHaloInFOFgroup);

        fprintf(fp,
                "%6"PRId64" %6d %6d  %6d %6d "
                "%6d %8d "
                "%8.3lf %8lld\n",
                i,
                forest[i].SnapNum, forest[i].Descendant, forest[i].FirstProgenitor, forest[i].NextProgenitor,
                forest[i].FirstHaloInFOFgroup, forest[i].NextHaloInFOFgroup,
                forest[i].Mvir, forest[i].MostBoundID);
            
        /* fprintf(stderr,"%12"PRId64" %12.5lf   %12d   %12"PRId64" %12"PRId64"  %12"PRId64" %14.6lf " */
        /* /\*         /\\* "%14.6lf %14.6lf %12.5lf %12.5lf %12.5lf " *\\/ *\/ */
        /* /\*         /\\* "%12.4lf %12.4lf %12.4lf " *\\/ *\/ */
        /* /\*         /\\* "%12.4lf %12.4lf %12.4lf " *\\/ *\/ */
        /* /\*         /\\* "%12.4lf %12.4lf %12d " *\\/ *\/ */
        /* /\*         "%12d %12d %12d \n", *\/ */
        /*         "%12d %12d \n", */
        /*         i, */
        /*         info[i].scale, forest[i].SnapNum, info[i].id, info[i].pid, info[i].upid, forest[i].Mvir, */
        /* /\*         /\\* forest[i].VelDisp, forest[i].Vmax, forest[i].Pos[0], forest[i].Pos[1], forest[i].Pos[2], *\\/ *\/ */
        /* /\*         /\\* forest[i].Vel[0], forest[i].Vel[1], forest[i].Vel[2], *\\/ *\/ */
        /* /\*         /\\* forest[i].Spin[0], forest[i].Spin[1], forest[i].Spin[2], *\\/ *\/ */
        /* /\*         /\\* forest[i].M_Mean200, forest[i].M_TopHat, forest[i].Len, *\\/ *\/ */
        /* /\*         forest[i].FirstProgenitor, forest[i].NextProgenitor, forest[i].Descendant, *\/ */
        /*         forest[i].FirstHaloInFOFgroup, forest[i].NextHaloInFOFgroup); */

        XASSERT( ( ! (forest[i].SnapNum == max_snapnum && forest[i].Descendant != -1)),
                 "Halo (forest[%"PRId64"] at last snapshot = %d has descendant = %d\n"
                 "Assigning descendants code needs to be checked\n",
                 i, max_snapnum, forest[i].Descendant);
                     
        if(forest[i].Descendant == -1) {
            //check that nothing points to it
            for(int64_t aa=0;aa<totnhalos;aa++) {
                XASSERT((forest[aa].FirstProgenitor != i) ||
                        (forest[aa].NextProgenitor  != i),
                        "i=%"PRId64" has no descendant but is assigned as firstprog or nextprog\n"
                        "forest[%"PRId64"].FirstProgenitor = %d "
                        "forest[%"PRId64"].NextProgenitor = %d \n",
                        i, aa, forest[aa].FirstProgenitor, aa, forest[aa].NextProgenitor);
                            
                            
            }
        } else {
            //has a descendant. make sure the correct halo points to it. either as a firstprog or
            //by following the chain of next prog.
            //then ensure that nothing else points to it.
            int expected_halo = forest[i].Descendant;
            //First the basic check that the descendant id corresponds to the one in info
            XASSERT(info[expected_halo].id == info[i].descid,
                    "Descendant ids must match\n"
                    "info[%d].id = %"PRId64" while info[%"PRId64"] = %"PRId64"\n",
                    expected_halo, info[expected_halo].id, i, info[i].descid);
                
            if(forest[expected_halo].FirstProgenitor != i) {
                //find the halo that points to `i'
                expected_halo = forest[expected_halo].FirstProgenitor;
                while(forest[expected_halo].NextProgenitor != -1 && forest[expected_halo].NextProgenitor != i) {
                    expected_halo = forest[expected_halo].NextProgenitor;
                    XASSERT(expected_halo >=0 && expected_halo < totnhalos,
                            "Invalid location for the source halo that points to %"PRId64"\n"
                            "Descendant location = %d. expected_halo = %d\n",
                            i, forest[i].Descendant, expected_halo);
                }
                XASSERT(forest[expected_halo].NextProgenitor == i,
                        "Should have found the source halo that points to %"PRId64"\n"
                        "Descendant location = %d\n"
                        "forest[%d].NextProgenitor, FirstProg  = %d,%d \n",
                        i, forest[i].Descendant, expected_halo,
                        forest[expected_halo].NextProgenitor,
                        forest[expected_halo].FirstProgenitor);
            } else {
                expected_halo = forest[i].Descendant;//Just to be explicit. This line will *NOT* change the value
            }

            for(int64_t aa=0;aa<totnhalos;aa++) {
                if(aa == expected_halo) {
                    continue;
                }

                //NO other halo should point to i (loop flow should not reach here for expected_halo)
                XASSERT( !  ((forest[aa].FirstProgenitor == i) ||
                             (forest[aa].NextProgenitor == i)),
                         "i=%"PRId64" should only be pointed to by %d\n"
                         "BUT: forest[%"PRId64"].FirstProgenitor = %d "
                         "forest[%"PRId64"].NextProgenitor = %d \n",
                         i, expected_halo, aa, forest[aa].FirstProgenitor, aa, forest[aa].NextProgenitor);
            }
        }

        /* ok now check the firsthaloinfofgroup + nexthaloinfofgroup*/
        if(info[i].pid !=-1) {
            continue;
        }

        //Only FOF halos will come here
        {
            int64_t aa = i+1;
            const int64_t fof_id = info[i].id;
            int64_t max_aa_inc = i;
            while(aa < totnhalos && info[aa].pid != -1) {
                XASSERT(forest[aa].FirstHaloInFOFgroup == i,
                        "forest[%"PRId64"].FirstHaloInFOFgroup=%d should point to %"PRId64" as fof halo\n",
                        aa, forest[aa].FirstHaloInFOFgroup, i);
                XASSERT(info[aa].upid == fof_id,
                        "info[%"PRId64"].upid = %"PRId64" instead of %"PRId64"\n",
                        aa, info[aa].upid, fof_id);
                if(aa > max_aa_inc) {
                    max_aa_inc = aa;
                }
                aa++;
            }

            aa = i;
            int64_t max_aa_next = i;
            while(aa >= 0 && aa < totnhalos) {
                XASSERT(forest[aa].FirstHaloInFOFgroup == i,
                        "forest[%"PRId64"].FirstHaloInFOFgroup=%d should point to %"PRId64" as fof halo\n",
                        aa, forest[aa].FirstHaloInFOFgroup, i);
                aa = forest[aa].NextHaloInFOFgroup;
                if(aa > max_aa_next) {
                    max_aa_next = aa;
                }
            }

            XASSERT(max_aa_next == max_aa_inc,
                    "All (sub) halos should consecutive\n"
                    "forest[%"PRId64"].FirstHaloInFOFgroup, Next = %d,%d \n"
                    "max located via next = %"PRId64" max located by increment = %"PRId64"\n",
                    i,
                    forest[i].FirstHaloInFOFgroup, forest[i].NextHaloInFOFgroup,
                    max_aa_next, max_aa_inc);

            for(int64_t bb=0;bb<totnhalos;bb++) {
                if(bb >= i || bb < max_aa_next) {
                    continue;
                }

                XASSERT(forest[bb].FirstHaloInFOFgroup != i,
                        "NO other halo should claim %"PRId64" as the FOF halo\n",
                        bb);
                XASSERT(forest[bb].NextHaloInFOFgroup != i,
                        "NO other halo should point to %"PRId64" as the next halo in FOF group\n",
                        bb);
            }
        }
    }

    fclose(fp);
}


int run_checks_on_new_locations(const struct locations *new_locations, const struct locations *locations, const int64_t ntrees)
{
    for(int64_t i=0;i<ntrees;i++){
        XPRINT(locations[i].forestid == new_locations[i].forestid,
               ANSI_COLOR_RED"locations[%"PRId64"].forestid = %"PRId64" does not equal new_locations[%"PRId64"].forestid = %"PRId64 ANSI_COLOR_RESET"\n",
               i, locations[i].forestid, i, new_locations[i].forestid);
        XPRINT(locations[i].tree_root == new_locations[i].tree_root,
               ANSI_COLOR_RED"locations[%"PRId64"].tree_root = %"PRId64" does not equal new_locations[%"PRId64"].tree_root = %"PRId64 ANSI_COLOR_RESET"\n",
               i, locations[i].tree_root, i, new_locations[i].tree_root);
        XPRINT(locations[i].fileid == new_locations[i].fileid,
               ANSI_COLOR_RED"locations[%"PRId64"].fileid = %"PRId64" does not equal new_locations[%"PRId64"].fileid = %"PRId64 ANSI_COLOR_RESET"\n",
               i, locations[i].fileid, i, new_locations[i].fileid);
        XPRINT(locations[i].offset == new_locations[i].offset,
               ANSI_COLOR_RED"locations[%"PRId64"].offset = %"PRId64" does not equal new_locations[%"PRId64"].offset = %"PRId64 ANSI_COLOR_RESET"\n",
               i, locations[i].offset, i, new_locations[i].offset);
               
        if(locations[i].forestid  != new_locations[i].forestid) return EXIT_FAILURE;
        if(locations[i].tree_root != new_locations[i].tree_root) return EXIT_FAILURE;
        if(locations[i].fileid    != new_locations[i].fileid) return EXIT_FAILURE;
        if(locations[i].offset    != new_locations[i].offset) return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

    
