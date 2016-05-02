#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "convert_trees_to_lhalo.h"
#include "progressbar.h"

void usage(int argc, char **argv)
{
    (void) argc;
    fprintf(stderr,"USAGE: %s <input consistent-trees directory> <output LHALOTREE directory> <particle mass (10^10 Msun/h units) >\n",
            argv[0]);
}    

int main(int argc, char **argv)
{
    char *input_dir, *output_dir;
    double part_mass = 0.0;
    char buffer[MAXLEN];
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
    const float inv_part_mass = 1.0f/part_mass;

    {
        const size_t expected_struct_size = 104;
        XASSERT(sizeof(struct output_dtype) == expected_struct_size,
                "sizeof output struct must exactly equal %zu bytes\n",
                expected_struct_size);
    }

    struct rlimit rlp;
    getrlimit(RLIMIT_NOFILE, &rlp);
    rlp.rlim_cur = rlp.rlim_max;
    setrlimit(RLIMIT_NOFILE, &rlp);
    
    struct timeval tstart, tend, t0, t1;
    gettimeofday(&tstart, NULL);
    char locations_filename[MAXLEN], forests_filename[MAXLEN];
    int64_t *forests=NULL, *tree_roots=NULL;
    my_snprintf(locations_filename, MAXLEN, "%s/locations.dat", input_dir);
    my_snprintf(forests_filename, MAXLEN, "%s/forests.list", input_dir);
    gettimeofday(&t0,NULL);
    fprintf(stderr, ANSI_COLOR_MAGENTA"Reading forests...."ANSI_COLOR_RESET"\n");
    const int64_t ntrees = read_forests(forests_filename, output_dir, &forests, &tree_roots);
    gettimeofday(&t1,NULL);
    fprintf(stderr, ANSI_COLOR_GREEN"Reading forests......done. Ntrees = %"PRId64". Time = %12.3lf seconds"ANSI_COLOR_RESET"\n\n", ntrees, ADD_DIFF_TIME(t0, t1));

    struct locations *locations = my_malloc(sizeof(*locations), ntrees);
    int64_t nfiles = 0, BOX_DIVISIONS=0;
    gettimeofday(&t0,NULL);
    fprintf(stderr, ANSI_COLOR_MAGENTA"Reading locations...."ANSI_COLOR_RESET"\n");
    const int64_t ntrees_loc = read_locations(locations_filename, ntrees, locations, &nfiles, &BOX_DIVISIONS);
    gettimeofday(&t1,NULL);
    fprintf(stderr, ANSI_COLOR_GREEN"Reading locations......done. Time = %12.3lf seconds"ANSI_COLOR_RESET"\n\n", ADD_DIFF_TIME(t0, t1));
    XASSERT(ntrees == ntrees_loc,
            "ntrees=%"PRId64" should be equal to ntrees_loc=%"PRId64"\n",
            ntrees, ntrees_loc);    

    /* the following function will sort locations and forests based on tree root id*/
    fprintf(stderr, ANSI_COLOR_MAGENTA"Assigning forest ids...."ANSI_COLOR_RESET"\n");
    gettimeofday(&t0,NULL);
    assign_forest_ids(ntrees, locations, forests, tree_roots);
    gettimeofday(&t1,NULL);
    fprintf(stderr, ANSI_COLOR_GREEN"Assigning forest ids.......done. Time = %12.3lf seconds"ANSI_COLOR_RESET"\n\n", ADD_DIFF_TIME(t0, t1));
    
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
    for (int i=0; i<BOX_DIVISIONS; i++) {
        for (int j=0; j<BOX_DIVISIONS; j++) {
            for(int k=0; k<BOX_DIVISIONS; k++) {
                my_snprintf(buffer,MAXLEN,"%s/tree_%d_%d_%d.dat", input_dir, i, j, k);
                int64_t id  = i*BOX_DIVISIONS*BOX_DIVISIONS + j*BOX_DIVISIONS + k;
                tree_inputs[id]  = my_fopen(buffer, "r");
                /* assert(setvbuf(tree_inputs[id], NULL, _IONBF, 0) == 0); */
                my_fseek(tree_inputs[id],0L, SEEK_END);
                inp_file_sizes[id] = ftello(tree_inputs[id]);
                rewind(tree_inputs[id]);

                tree_inputs_fd[id]  = fileno(tree_inputs[id]);

                my_snprintf(buffer,MAXLEN,"%s/lhalotree.bin.%"PRId64"", output_dir, id);
                unlink(buffer);
                my_snprintf(tree_outputs_fname[id], LOCATIONS_FILENAME_SIZE, "lhalotree.bin.%"PRId64"",id);
                tree_outputs[id] = my_fopen(buffer, "w");
                tree_outputs_fd[id] = fileno(tree_outputs[id]);
                tree_outputs_fd_offset[id] = 0;
            }
        }
    }


    /* the following function will sort locations based on 1) filename 2) offsets */
    fprintf(stderr, ANSI_COLOR_MAGENTA"Sorting locations based on file offsets...."ANSI_COLOR_RESET"\n");
    sort_locations_file_offset(ntrees, locations);
    fprintf(stderr, ANSI_COLOR_GREEN"Sorting locations based on file offsets........done"ANSI_COLOR_RESET"\n\n");

    int interrupted=0;//for progressbar
 
    /* holder to check later that bytes have been assigned */
    for(int64_t i=0;i<ntrees;i++) {
        locations[i].bytes = -1;/* Make sure bytes is a signed type! */
    }
    /* Create a copy of current locations */    
    struct locations *output_locations;
    char locations_bin_filename[MAXLEN];
    my_snprintf(locations_bin_filename,MAXLEN,"%s/locations.bin", output_dir);    
    FILE *locations_binary_fp = fopen(locations_bin_filename,"r");
    int compute_bytes = 1;

    //There is a locations file but may be it was from a different run or different data-set
    //This will be wasteful but allocate a new locations and then read in the data and check
    //that everything is okay. Only then, trust the number of bytes. 
    if(locations_binary_fp != NULL) {
        // Found the locations.bin file -> read it in and avoid computing the number of bytes
        gettimeofday(&t0, NULL);
        fprintf(stderr,ANSI_COLOR_MAGENTA"Reading binary locations file "ANSI_COLOR_GREEN"`%s'"ANSI_COLOR_MAGENTA"..."ANSI_COLOR_RESET"\n",locations_bin_filename);
        size_t dummy;
        int64_t ntrees_in_file;
        struct locations *tmp_locations = my_malloc(sizeof(*tmp_locations), ntrees);
        my_fread(&dummy, sizeof(dummy), 1, locations_binary_fp);
        assert(dummy == sizeof(*tmp_locations));
        my_fread(&ntrees_in_file, sizeof(ntrees_in_file), 1, locations_binary_fp);
        assert(ntrees_in_file == ntrees);
        my_fread(tmp_locations, sizeof(*tmp_locations), ntrees_in_file, locations_binary_fp);
        fclose(locations_binary_fp);
        locations_binary_fp = NULL;

        //If all the fields that should agree, (essentially everything set in locations so far),
        //then we can be sure we are not corrupting the data. 
        const int status = run_checks_on_new_locations(tmp_locations, locations, ntrees);
        if (status == EXIT_SUCCESS) {
            output_locations = my_malloc(sizeof(*output_locations), ntrees);
            assert(sizeof(*output_locations) == sizeof(*tmp_locations) && "locations struct is varying in size! The sky is falling!!");
            memcpy(output_locations, tmp_locations, sizeof(*tmp_locations) * ntrees);
            gettimeofday(&t1, NULL);
            fprintf(stderr,ANSI_COLOR_MAGENTA"Reading binary locations file "ANSI_COLOR_GREEN"`%s'"ANSI_COLOR_MAGENTA".....done. Time = %12.3lf seconds"ANSI_COLOR_RESET"\n\n",
                    locations_bin_filename, ADD_DIFF_TIME(t0, t1));
            compute_bytes = 0;
            fprintf(stderr,ANSI_COLOR_MAGENTA"If you see parse errors later on, that could be because the `forests.bin' and `locations.bin' do not correspond to `%s' "ANSI_COLOR_RESET"\n", input_dir);
            fprintf(stderr,ANSI_COLOR_MAGENTA"In that case, delete the `forests.bin' and `locations.bin' files in `%s' and restart the code"ANSI_COLOR_RESET"\n",output_dir);
        } else {
            fprintf(stderr,ANSI_COLOR_RED"ERROR: Locations in file `%s' does not agree with currently read in locations. Deleting the `locations.bin'.."ANSI_COLOR_RESET"\n",locations_bin_filename);
            unlink(locations_bin_filename);
            fprintf(stderr,ANSI_COLOR_RED"Forests data might have been corrupted as well. Deleting the `forests.bin' file in output directory"ANSI_COLOR_RESET"\n");
            char forest_bin_file[MAXLEN];
            my_snprintf(forest_bin_file, MAXLEN, "%s/forests.bin", output_dir);
            unlink(forest_bin_file);

            //Ideally, I should restart the code but that will require goto's.
            //Just exit and tell the user to restart the code. 
            fprintf(stderr,ANSI_COLOR_RED"\tPlease ensure the `%s/forests.bin' and `%s/locations.bin' files have indeed been deleted"ANSI_COLOR_RESET"\n", output_dir, output_dir);
            fprintf(stderr,ANSI_COLOR_GREEN"\tAfter that, please restart the code with `%s %s %s %s'"ANSI_COLOR_RESET"\n",argv[0], argv[1], argv[2], argv[3]);
            fprintf(stderr,"Exiting now...\n");
            exit(EXIT_FAILURE);
        }
    }

        
    //Is there a previous run that I could read in?
    if(compute_bytes == 1) {
        output_locations = my_malloc(sizeof(*output_locations), ntrees);
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
        interrupted=0;
        init_my_progressbar(ntrees, &interrupted);
        
        for(int64_t i=1;i<=ntrees-1;i++) {
            my_progressbar(i, &interrupted);
            const int64_t fileid = locations[i].fileid;
            
            /* Are we starting on a new file ?*/
            if(start_fileid != fileid) {
                /* fill out the bytes for the last tree in the previous file */
                const int64_t num_bytes = compute_numbytes_with_off(inp_file_sizes[start_fileid], start);
                locations[i-1].bytes = num_bytes;
                output_locations[i-1].bytes = num_bytes;
                
                /* now we reset the start fields */
                start = locations[i].offset;
                start_fileid = locations[i].fileid;
                continue;
            }
            const int64_t current_offset_guess = locations[i].offset - guess_max_linesize;
            
#if 1
            my_fseek(tree_inputs[fileid], current_offset_guess, SEEK_SET);
            while(1) {
                const int a = fgetc(tree_inputs[fileid]);
                if(a == EOF) {
                    fprintf(stderr,"Encountered EOF while looking for end of current tree\n");
                    exit(EXIT_FAILURE);
                }
                const unsigned char c = (unsigned char) a;
                if(c == '\n') {
                    //Why is this start_fileid rather than fileid?
                    const int64_t num_bytes = compute_numbytes(tree_inputs[start_fileid], start);
                    locations[i-1].bytes = num_bytes;
                    output_locations[i-1].bytes = num_bytes;
                    /* fprintf(stderr,"%"PRId64"\n",num_bytes); */
                    start = locations[i].offset;
                    break;
                }
            }
            
#else
            assert(MAXLEN > guess_max_linesize);
            int64_t curr_offset = current_offset_guess;
            int64_t bytes_this_read = (int64_t) pread(tree_inputs_fd[start_fileid], buffer, guess_max_linesize, current_offset_guess);
            assert(bytes_this_read == guess_max_linesize);
            {
                int found = 0;
                for(int64_t ii=0;ii<bytes_this_read; ii++) {
                    curr_offset++;
                    if(buffer[ii] == '\n'){
                        found = 1;
                        const int64_t num_bytes = compute_numbytes_with_off(curr_offset, start);
                        locations[i-1].bytes = num_bytes;
                        output_locations[i-1].bytes = num_bytes;
                        start = locations[i].offset;
                        break;
                    } 
                }
                assert(found == 1);
            }
#endif
        }
    

    
        /* fill out the bytes for the last tree */
        {
            start = locations[ntrees-1].offset;
            const int64_t fileid = locations[ntrees-1].fileid;
            
#if 0        
            my_fseek(tree_inputs[fileid], 0L, SEEK_END);
            const int64_t num_bytes = compute_numbytes(tree_inputs[fileid], start);
#else
            const int64_t num_bytes = compute_numbytes_with_off(inp_file_sizes[fileid], start);
#endif        
            
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
            
            for(int64_t i=0;i<nfiles;i++) {
                XASSERT(total_tree_bytes[i] < inp_file_sizes[i],
                        "Bytes in tree = %"PRId64" must be smaller than file size = %"PRId64"\n",
                        total_tree_bytes[i], inp_file_sizes[i]);
            }
            free(total_tree_bytes);
            
            /* Output this locations data */
            assert(locations_binary_fp == NULL);
            locations_binary_fp = my_fopen(locations_bin_filename, "w");
            size_t dummy = sizeof(*locations);
            my_fwrite(&dummy, sizeof(dummy), 1, locations_binary_fp);
            my_fwrite(&ntrees, sizeof(ntrees), 1, locations_binary_fp);
            my_fwrite(locations, sizeof(*locations), ntrees, locations_binary_fp);//totnforests in this file 
            fclose(locations_binary_fp);

        }
    } //computing the number of bytes and saving the info as locations.bin file
    

    
    /* Now assign all trees in the same forest to the same file. 
       The output fileids goes into output_locations (which is otherwise a copy of locations).

       Both locations and output_locations are sorted by ForestID, FileID, Offset (in that order).
     */
    struct forest_info *forest_info = assign_trees_in_forest_to_same_file(ntrees, locations, output_locations, nfiles);
    int64_t nforests = forest_info->nforests;//should really be a const but I am changing nforests later for debugging. 
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

/* #define RESET_FORESTS  (3) */
    
    for(int64_t i=0;i<nfiles;i++) {
        FILE *fp = tree_outputs[i];
        rewind(fp);
        const int zero = 0;
        my_fwrite(&zero, sizeof(int), 1, fp);//totnforests in this file 
        my_fwrite(&zero, sizeof(int), 1, fp);//totnhalos in this file
#ifdef RESET_FORESTS        
        totnforests_per_file[i] = RESET_FORESTS;
#endif        
        for(int64_t j=0;j<totnforests_per_file[i];j++) {
            my_fwrite(&zero, sizeof(int), 1, fp);//one zero for each forest in this file. will contain the number of halos in forest.
        }
        fflush(fp);
        fsync(fileno(fp));/* ensure data really are flushed to disk */
        offsets_per_file[i] = ftello(fp);
    }
    
#ifdef RESET_FORESTS
    nforests = RESET_FORESTS;
#endif    

    fprintf(stderr, ANSI_COLOR_MAGENTA"Writing out (nforests=%"PRId64") in LHALOTREE format...."ANSI_COLOR_RESET"\n", nforests);
    interrupted=0;
    init_my_progressbar(nforests, &interrupted);

    /* Now copy each one of the forests. Note a forest can have multiple trees */
    int64_t nhalos_allocated = 1000000;//1 million halos
    struct output_dtype *forest_halos  = my_malloc(sizeof(*forest_halos), nhalos_allocated);
    struct additional_info *info = my_malloc(sizeof(*info), nhalos_allocated);
    int64_t tree_index = 0;
    my_snprintf(buffer, MAXLEN, "%s/output_order_forests.dat", output_dir);
    FILE *output_order_locations = my_fopen(buffer, "w");
    fprintf(output_order_locations,"#TreeRootID ForestID\n");

    int64_t total_num_halos_written_all_files = 0;
    for(int64_t i=0;i<nforests;i++) {
        my_progressbar(i, &interrupted);
        int64_t forest_offset = 0;
        int64_t totnhalos = 0;
        /* fprintf(stdout, "%"PRId64" %"PRId64" ", i, forest_info->num_trees[i]); */
        for(int64_t j=0;j<forest_info->num_trees[i];j++) {
            fprintf(output_order_locations,"%"PRId64" %"PRId64"\n",
                    locations[tree_index].tree_root, locations[tree_index].forestid);
            fflush(output_order_locations);
            int64_t fileid = locations[tree_index].fileid;
            const int64_t nhalos = read_tree_into_forest(&nhalos_allocated, &forest_halos, forest_offset, &info,
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
            totnhalos += nhalos;
            total_num_halos_written_all_files += nhalos;
        }
        
        const int64_t out_fileid = forest_info->fileid[i];
        const int forestindex_thisfile = nforests_written_per_file[out_fileid];
        nhalos_per_forest_per_file[out_fileid][forestindex_thisfile] = totnhalos;
        totnhalos_per_file[out_fileid] += totnhalos;
        nforests_written_per_file[out_fileid]++;
        nhalos_written_per_file[out_fileid] += totnhalos;
        /* const int verbose = (i == 8780) ? 1:0; */
        int verbose = 0;

        /* Fix flybys -> multiple roots at z=0 must be joined such that only one root remains */
        int status = fix_flybys(totnhalos, forest_halos, info, verbose);
        if(status != EXIT_SUCCESS) {
            fprintf(stderr,ANSI_COLOR_RED"ERROR while trying to convert Forest id = %"PRId64" with ntrees = %"PRId64 " last tree index = %10"PRId64 ANSI_COLOR_RESET"\n",
                    forest_info->forestid[i], forest_info->num_trees[i], tree_index);

            tree_index -= forest_info->num_trees[i];
            for(int64_t j=0;j<forest_info->num_trees[i];j++) {
                fprintf(stderr,"Tree Root = %10"PRId64" inp. fileid = %4"PRId64" output fileid = %4"PRId64"\n",
                        locations[tree_index].tree_root,
                        locations[tree_index].fileid,
                        output_locations[tree_index].fileid);
                XASSERT(locations[tree_index].forestid == output_locations[tree_index].forestid,
                        "locations[%"PRId64"].forestid = %"PRId64" must equal output_locations[%"PRId64"].forestid = %"PRId64"\n",
                        tree_index,locations[tree_index].forestid, tree_index, output_locations[tree_index].forestid);
                tree_index++;
            }        
            exit(EXIT_FAILURE);
                    
        }

        /* Entire tree is loaded in. Fix upid's*/
        const int max_snapnum = fix_upid(totnhalos, forest_halos, info, &interrupted, verbose);
        
        /* Now the entire tree is loaded in. Assign the mergertree indices */
        assign_mergertree_indices(totnhalos, forest_halos, info, max_snapnum);

        const int64_t num_bytes = sizeof(struct output_dtype) * totnhalos;
        forest_info->num_binary_bytes[i]  = num_bytes;
        forest_info->offset[i]    = offsets_per_file[out_fileid];
        forest_info->num_halos[i]  = totnhalos;
        offsets_per_file[out_fileid] += num_bytes;
        
#ifdef RESET_FORESTS
        validate_fields(totnhalos, forest_halos, info, max_snapnum);
#endif
        /* write out the forest in binary */
        my_fwrite(forest_halos, sizeof(*forest_halos), totnhalos, tree_outputs[out_fileid]);
        
    }
    fclose(output_order_locations);

    for(int64_t i=0;i<nfiles;i++) {
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
                fprintf(stderr,"\nIn file `%s' i = %"PRId64" nforests_written = %d does not agree with (int64_t) totnforests_per_file = %"PRId64"\n"
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
                fprintf(stderr,"\nIn file `%s' i = %"PRId64" nhalos_written = %d does not agree with (int64_t) totnhalos_per_file = %"PRId64"\n"
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
        fprintf(fp,"#######################################################################################################################################\n");
        fprintf(fp,"#         ForestID     FileID         Filename                  Offset            BinaryBytes         AsciiBytes            NumHalos   \n");
        fprintf(fp,"#######################################################################################################################################\n");

        for(int64_t i=0;i<forest_info->nforests;i++) {
            fprintf(fp,"%18"PRId64" %10"PRId64"      %s   %18"PRId64"  %18"PRId64"  %18"PRId64"    %18"PRId64"\n",
                    forest_info->forestid[i],forest_info->fileid[i], forest_info->filename[i],
                    forest_info->offset[i], forest_info->num_binary_bytes[i], forest_info->num_ascii_bytes[i],forest_info->num_halos[i]);
        }
        fclose(fp);
    }

    
    /* close open file pointers + free memory for file pointers */
    for(int64_t i=0;i<nfiles;i++) {
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
    free(forest_info->num_ascii_bytes);
    free(forest_info->num_binary_bytes);
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
    free(forest_halos);free(info);
    
    gettimeofday(&tend, NULL);
    fprintf(stderr,"\n\nWrote out %"PRId64" halos (across all files) contained in %"PRId64" trees in %"PRId64" forests in the LHALOTree format. Time taken = %0.2g seconds\n\n",
            total_num_halos_written_all_files, ntrees, nforests, ADD_DIFF_TIME(tstart, tend));
    
    return EXIT_SUCCESS;
}
    
    
