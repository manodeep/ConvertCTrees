#include <stdio.h>
#include <stdlib.h>

#include "../output_datatype.h"
#include "../sglib.h"
#include "../utils.h"

#define _FILE_OFFSET_BITS 64

void usage(int argc, char **argv)
{
    (void) argc;
    fprintf(stderr,"USAGE: `%s' <first lhalotree file> <second lhalotree file> <check mergertree indices>\n",argv[0]);
}    

struct output_dtype * read_lhalotree(const char *filename, int *ntrees, int *totnhalos, int **nhalos_per_tree)
{
    FILE *fp = my_fopen(filename,"r");
    my_fread(ntrees, sizeof(*ntrees), 1, fp);
    my_fread(totnhalos, sizeof(*totnhalos), 1, fp);
    *nhalos_per_tree = my_malloc(sizeof(**nhalos_per_tree), *ntrees);
    my_fread(*nhalos_per_tree, sizeof(**nhalos_per_tree), *ntrees, fp);
    struct output_dtype *all_trees = my_malloc(sizeof(*all_trees), *totnhalos);
    my_fread(all_trees, sizeof(*all_trees), *totnhalos, fp);

    const off_t curpos = ftello(fp);
    my_fseek(fp, 0, SEEK_END);
    const off_t endpos = ftello(fp);
    XASSERT(endpos == curpos,
            "Managed to read in entire file."
            "Position after reading in all data  = %"PRId64" FILE endpos = %"PRId64"\n",
            curpos, endpos);
    fclose(fp);

    return all_trees;
}    

int main(int argc, char **argv)
{
    if(argc != 4) {
        usage(argc, argv);
        fprintf(stderr,"exiting\n");
        exit(EXIT_FAILURE);
    }

    char *file1 = argv[1];
    char *file2 = argv[2];
    const int check_mergertree_indices = atoi(argv[3]);
    int first_ntrees=0, first_nhalos = 0;
    int second_ntrees=0, second_nhalos=0;
    int *first_nhalos_per_tree=NULL, *second_nhalos_per_tree=NULL;
    struct output_dtype *first = read_lhalotree(file1, &first_ntrees, &first_nhalos, &first_nhalos_per_tree);
    struct output_dtype *second = read_lhalotree(file2, &second_ntrees, &second_nhalos, &second_nhalos_per_tree);
    const int LEN_EPS = 2;
    const float MASS_EPS = 1e-3, POS_EPS = 1e-3;
    
    XASSERT(first_ntrees == second_ntrees,
            "Two files must have same number of trees\n"
            "File `%s' contains  : %d trees\n"
            "File `%s' contains  : %d trees\n",
            file1, first_ntrees,
            file2, second_ntrees);

    XASSERT(first_nhalos == second_nhalos,
            "Two files must have same number of halos\n"
            "File `%s' contains  : %d halos\n"
            "File `%s' contains  : %d halos\n",
            file1, first_nhalos,
            file2, second_nhalos);

    /* If control flow reaches here, then the total number of halos and number of trees must be identical
       Save some typing :) 
     */
    const int ntrees = first_ntrees;
    const int totnhalos = first_nhalos;
    for(int i=0;i<ntrees;i++) {
        XASSERT(first_nhalos_per_tree[i] == second_nhalos_per_tree[i],
                "i=%d nhalos1=%d nhalos2=%d\n",
                i,first_nhalos_per_tree[i],second_nhalos_per_tree[i]);
    }

    
    /* Okay so even the number of halos per tree are good. Now check individual halos.
       Sort them by (Ctrees) HaloID (which has been conveniently stored in MostBoundID
       by the code in the parent folder. 
     */

    //Sort the first list of halos
    int max_snapnum = -1;
    for(int i=0;i<totnhalos;i++) {
        if(first[i].SnapNum > max_snapnum) {
            max_snapnum = first[i].SnapNum;
        }
    }
    
    int offset = 0;
    for(int itree=0;itree<ntrees;itree++) {
        const int nhalos = first_nhalos_per_tree[itree];
        struct output_dtype *first_halos  = &(first[offset]);
        struct output_dtype *second_halos = &(second[offset]);

        if(check_mergertree_indices == 0) {
#define MOSTBOUND_COMPARATOR(x, y)   (((x.MostBoundID)>(y.MostBoundID)?1:((x.MostBoundID)<(y.MostBoundID)?-1:0)))
            SGLIB_ARRAY_SINGLE_HEAP_SORT(struct output_dtype, first_halos, nhalos, MOSTBOUND_COMPARATOR);
            SGLIB_ARRAY_SINGLE_HEAP_SORT(struct output_dtype, second_halos, nhalos, MOSTBOUND_COMPARATOR);
#undef MOSTBOUND_COMPARATOR            
        }
        
        for(int i=0;i<nhalos;i++) {
            XPRINT(first_halos->MostBoundID == second_halos->MostBoundID,
                   "itree = %d, i = %d\n"
                   "first->id  = %lld (`%s')\n"
                   "second->id = %lld (`%s')\n",
                   itree, i,
                   first_halos->MostBoundID, file1,
                   second_halos->MostBoundID, file2);
            
            XPRINT(abs(first_halos->Len - second_halos->Len) <= LEN_EPS,
                   "itree = %d, i = %d\n"
                   "first->Len  = %d (`%s')\n"
                   "second->Len = %d (`%s')\n",
                   itree, i,
                   first_halos->Len, file1,
                   second_halos->Len, file2);
            
            XPRINT(fabsf(first_halos->Mvir - second_halos->Mvir) <= MASS_EPS,
                   "itree = %d, i = %d\n"
                   "first->Mvir  = %lf (`%s')\n"
                   "second->Mvir = %lf (`%s')\n",
                   itree, i,
                   first_halos->Mvir, file1,
                   second_halos->Mvir, file2);
            
            XPRINT(fabsf(first_halos->M_Mean200 - second_halos->M_Mean200) <= MASS_EPS,
                    "itree = %d, i = %d\n"
                    "halo id = %lld snapnum = %d\n"
                    "first->M_Mean200  = %lf (`%s')\n"
                    "second->M_Mean200 = %lf (`%s')\n",
                    itree, i,
                    first_halos->MostBoundID, first_halos->SnapNum,
                    first_halos->M_Mean200, file1,
                    second_halos->M_Mean200, file2);
            
            XPRINT(fabsf(first_halos->M_TopHat - second_halos->M_TopHat) <= MASS_EPS,
                    "itree = %d, i = %d\n"
                    "halo id = %lld snapnum = %d\n"
                    "first->M_TopHat  = %lf (`%s')\n"
                    "second->M_TopHat = %lf (`%s')\n",
                    itree, i,
                    first_halos->MostBoundID, first_halos->SnapNum,
                    first_halos->M_TopHat, file1,
                    second_halos->M_TopHat, file2);
            
            for(int j=0;j<3;j++) {
                XPRINT(fabsf(first_halos->Pos[j] - second_halos->Pos[j]) <= POS_EPS,
                        "itree = %d, i = %d\n"
                        "first->Pos[j]  = %lf (`%s')\n"
                        "second->Pos[j] = %lf (`%s')\n",
                        itree, i,
                        first_halos->Pos[j], file1,
                        second_halos->Pos[j], file2);
                
                XPRINT(fabsf(first_halos->Vel[j] - second_halos->Vel[j]) <= POS_EPS,
                        "itree = %d, i = %d\n"
                        "first->Vel[j]  = %lf (`%s')\n"
                        "second->Vel[j] = %lf (`%s')\n",
                        itree, i,
                        first_halos->Vel[j], file1,
                        second_halos->Vel[j], file2);
                
                XPRINT(fabsf(first_halos->Spin[j] - second_halos->Spin[j]) <= POS_EPS,
                        "itree = %d, i = %d\n"
                        "first->Spin[j]  = %lf (`%s')\n"
                        "second->Spin[j] = %lf (`%s')\n",
                        itree, i,
                        first_halos->Spin[j], file1,
                        second_halos->Spin[j], file2);
            }
            
            XPRINT(fabsf(first_halos->VelDisp - second_halos->VelDisp) <= POS_EPS,
                    "itree = %d, i = %d\n"
                    "first->VelDisp  = %lf (`%s')\n"
                    "second->VelDisp = %lf (`%s')\n",
                    itree, i,
                    first_halos->VelDisp, file1,
                    second_halos->VelDisp, file2);
            
            XPRINT(fabsf(first_halos->Vmax - second_halos->Vmax) <= POS_EPS,
                    "itree = %d, i = %d\n"
                    "first->Vmax  = %lf (`%s')\n"
                    "second->Vmax = %lf (`%s')\n",
                    itree, i,
                    first_halos->Vmax, file1,
                    second_halos->Vmax, file2);
            
            XPRINT(first_halos->SnapNum == second_halos->SnapNum,
                    "itree = %d, i = %d\n"
                    "first->SnapNum  = %d (`%s')\n"
                    "second->SnapNum = %d (`%s')\n",
                    itree, i,
                    first_halos->SnapNum, file1,
                    second_halos->SnapNum, file2);
            
            if((first_halos->SnapNum == max_snapnum && first_halos->Descendant != -1) ||
               (second_halos->SnapNum == max_snapnum && second_halos->Descendant != -1)) {
                fprintf(stderr,"itree = %d, i = %d max_snap = %d\n"
                        "first->Descendant,snapshot  = %d, %d  (`%s')\n"
                        "second->Descendant,snapshot = %d, %d  (`%s')\n",
                        itree, i, max_snapnum,
                        first_halos->Descendant, first->SnapNum, file1,
                        second_halos->Descendant, second->SnapNum, file2);
            }

            if(check_mergertree_indices == 1) {
                if(first_halos->Descendant != -1 && second_halos->Descendant != -1 ) {
                    XPRINT((first_halos->Descendant >= 0 && first_halos->Descendant < nhalos),
                            "first->desc = %d is not within range [0, %d)\n",
                            first_halos->Descendant, nhalos);
                    XPRINT((second_halos->Descendant >= 0 && second_halos->Descendant < nhalos),
                            "second->desc = %d is not within range [0, %d)\n",
                            second_halos->Descendant, nhalos);
                    
                    
                    XPRINT(first[offset + first_halos->Descendant].MostBoundID == second[offset + second_halos->Descendant].MostBoundID,
                           "itree = %d, i = %d nhalos = %d\n"
                           /* "first->Descendant, descid  = "ANSI_COLOR_RED"%d, %lld (`%s')"ANSI_COLOR_RESET"\n" */
                           /* "second->Descendant, descid = "ANSI_COLOR_RED"%d, %lld (`%s')"ANSI_COLOR_RESET"\n", */
                           "first->Descendant, descid, id  = %d, %lld, %lld (`%s')\n"
                           "second->Descendant, descid, id  = %d, %lld, %lld (`%s')\n",
                           itree, i, nhalos, 
                           first_halos->Descendant, first[offset + first_halos->Descendant].MostBoundID, first_halos->MostBoundID, file1,
                           second_halos->Descendant,  second[offset + second_halos->Descendant].MostBoundID, second_halos->MostBoundID, file2);
                }
                
                
                if(first_halos->FirstProgenitor != -1 && second_halos->FirstProgenitor != -1 ) {
                    XPRINT(first[offset + first_halos->FirstProgenitor].MostBoundID == second[offset + second_halos->FirstProgenitor].MostBoundID,
                           "itree = %d, i = %d\n"
                           "first->FirstProgenitor, firstprogid  = %d,%lld (`%s')\n"
                           "second->FirstProgenitor, firstprogid = %d,%lld (`%s')\n",
                           itree, i, 
                           first_halos->FirstProgenitor, first[offset + first_halos->FirstProgenitor].MostBoundID, file1,
                           second_halos->FirstProgenitor, second[offset + second_halos->FirstProgenitor].MostBoundID, file2);
                }
                
                
                if(first_halos->NextProgenitor != -1 && second_halos->NextProgenitor != -1 ) {
                    XPRINT(first[offset + first_halos->NextProgenitor].MostBoundID == second[offset + second_halos->NextProgenitor].MostBoundID,
                           "itree = %d, i = %d\n"
                           "first->NextProgenitor, nextprogid  = %d,%lld (`%s')\n"
                           "second->NextProgenitor, nextprogid = %d,%lld (`%s')\n",
                           itree, i, 
                           first_halos->NextProgenitor, first[offset + first_halos->NextProgenitor].MostBoundID, file1,
                           second_halos->NextProgenitor, second[offset + second_halos->NextProgenitor].MostBoundID, file2);
                }
                
                XPRINT(first_halos->FirstHaloInFOFgroup == second_halos->FirstHaloInFOFgroup,
                       "itree = %d, i = %d\n"
                       "first->FirstHaloInFOFgroup  = %d (`%s')\n"
                       "second->FirstHaloInFOFgroup = %d (`%s')\n",
                       itree, i,
                       first_halos->FirstHaloInFOFgroup, file1,
                       second_halos->FirstHaloInFOFgroup, file2);
                
                XPRINT(first_halos->NextHaloInFOFgroup == second_halos->NextHaloInFOFgroup,
                       "itree = %d, i = %d\n"
                       "first->NextHaloInFOFgroup  = %d (`%s')\n"
                       "second->NextHaloInFOFgroup = %d (`%s')\n",
                       itree, i,
                       first_halos->NextHaloInFOFgroup, file1,
                       second_halos->NextHaloInFOFgroup, file2);
            }
            
            first_halos++;
            second_halos++;
        }

        offset += nhalos;
    }
    free(first_nhalos_per_tree);
    free(second_nhalos_per_tree);
    
    free(first);
    free(second);
    return EXIT_SUCCESS;
}    
    
