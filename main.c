//
//  main.c
//  iGEMtoehold
//
//  Created by David Chu on 22/6/16.
//  Copyright © 2016 iGEM. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <sys/socket.h>

#include  "data_structures.h"
#include  "params.h"
#include  "utils.h"
#include  "eval.h"
#include  "fold.h"
#include  "part_func.h"

#define len  1001



void complement(char *string) {// make complement
    int i;
    for (i=0; i<strlen(string); i++) {
        if (string[i] == 'A')
            string[i] = 'U';
        else if (string[i] == 'C')
            string[i] = 'G';
        else if (string[i] == 'U')
            string[i] = 'A';
        else if (string[i] == 'G')
            string[i] = 'C';
    }
    
}

void makedna(char *string) {//U to T
    int i;
    for (i=0; i<strlen(string); i++) {
        if (string[i] == 'U')
            string[i] = 'T';
    }
    
}

void makerna(char *string) {//U to T
    int i;
    for (i=0; i<strlen(string); i++) {
        if (string[i] == 'T')
            string[i] = 'U';
    }
    
}

void Upcase(char str[1000]) {//unused
    int i;
    for (i = 0; i < 1000; i++) {
        //if (str[i] >= 'a' && str[i] <= 'z');
        str[i] = str[i] - 32;
    }
}


int check_valid(char str[1001]) {
    int i;
    int check = 0;
    for (i = 0; str[i] != '\0'; i++) {
        
        if (!(str[i] == 'A' || str[i] == 'U' || str[i] == 'C' || str[i] == 'G' || str[i] == '\n')) {
            //printf("%d \n",str[i]);
            check++;
        }
    }
    if (check>0){
        printf("Invalid Input :%d\n", check);
    }
    return check;
}


int choice_valid(int choice) {
    int check = 0;
    check = 0;
        if (choice<0 || choice>3) {
            //printf("%d \n",str[i]);
            check++;
        }
    if (check>0){
        printf("Invalid Input :%d\n", check);
    }
    return check;
}

void reversestring(char *string){//reverse the string
    int i, tmp;
    int j = strlen(string) - 1;
    for (i=0; i<j; i++) {
        tmp = string[i];
        string[i] = string[j];
        string[j] = tmp;
        j--;
        
    }
}

int look(char choice, int between, char arry[len], char rbs[500], int l, int window, int n, FILE *textfile) {
    char ecor1[] = "GAATTC"; //cutting site
    char not1[] = "GCGGCCGC"; //cutting site
    char xba1[] = "TCTAGA"; //cutting site
    char uaurbs[len] = "UAU"; //rbs site, starts with UAU
    char auarbs[len] = "AUA"; //rbs site, starts with AUA
    char link[] = "AACCUGGCGGCAGCGCAAAAG";//linker sequence + end cutting site
    char temp[31] = {'\0'};//store the 30bp RNA thing([31] to store the null character)
    char temp2[1000] = {'\0'};//store toehold switch sequence
    char temp3[100] = {'\0'};//store toehold switch DNA sequence
    char temp4[100] = {'\0'};//temporary storage for additional things
    char seq[100] = {'\0'};//temporary storage free energy process
    char dimer[10000] = {'\0'};//trigger switch dimer sequence
    char rbstemp[500] = {'\0'};;//turns all rbs things into '.'
    char rbsdimer[500] = {'\0'};;//rbs inside dimer
    char AUG[4] = "...";//chekc if AUG is "..."
    char AUGdimer[4] = "000";//AUG in dimer
    int j;
    int k;
    int x;
    int count = 0;// if count = 0, window++, successful switch makes count 1, unsuccessful makes 2;
    int count2 = 0;// to prevent window gain after consecutive unsuccessful switches
    int count3 = between;//space between each group
    int location = 0;// location
    int newwindow = 2;//new window or not
    int correctfolding = 0;//correct dimer folding rbs and aug
    double mfedifference = DBL_MIN;//difference between mfe for comparison
    double mfetemp = 0;//mfe of bestswitch
    char bestswitch[100] = {'\0'};//best switch with lowest mfe in a group
    int templocation = 0;//location of best switch
    
    //printf("started");
    
    strcpy(rbstemp, rbs);
    for (j=0; j<strlen(rbstemp); j++) {
        rbstemp[j] = '.';//turns all rbstemp things into '.'
    }
    
    strcat(uaurbs, rbs);
    strcat(uaurbs, "AUA");
    strcat(auarbs, rbs);
    strcat(auarbs, "UAU");
    
    if (l < 31){
        printf("Too short! ");
        return 0;
    }
    for (j = 0; j < (l - 31); j++) {
        if (count == 0){
            if (count2 == 0 && count3 == 3) {
                newwindow = 1;//newwindow true
                window++;//new window
                count3 = 0;
                count2 = 1;
                if((window>1) && (correctfolding == 1)){
                    correctfolding = 0;
                //start printing best switch of the group
                fprintf(textfile, "\n\nSwitch with highest mfe difference in group number %d\nLocation: %d\nSequence: ", (window-1), templocation);
                    for (k = 0; k < strlen(bestswitch); k++){//print final switch sequence in txt file
                        fprintf(textfile, "%c", bestswitch[k]);
                    }
                fprintf(textfile, "\nMFE: %fkcal/mol\nMFE difference: %fkcal/mol\n\n\n\n\n\n\n", mfetemp, mfedifference);
                //end printing best switch of the 
                }
            }
            else{
                newwindow = 0;//newwindow false
            }
            count = 2;
        }
        count = 0;//reset count
        count3++;
        location++;//location increase
        if (arry[j + 12] == 'A' || (arry[j + 12] == 'U' && arry[j + 13] != 'A' && arry[j + 13] != 'U' && arry[j + 14] != 'A' && arry[j + 14] != 'U')) {
            //printf("1's if! ");
            if (!((arry[j + 6] == 'U' && arry[j + 7] == 'A' && arry[j + 8] == 'A') || (arry[j + 6] == 'U' && arry[j + 7] == 'A' && arry[j + 8] == 'G')
                  || (arry[j + 6] == 'U' && arry[j + 7] == 'G' && arry[j + 8] == 'A') || (arry[j + 6] == 'U' && arry[j + 7] == 'C' && arry[j + 8] == 'A')
                  || (arry[j + 6] == 'U' && arry[j + 7] == 'U' && arry[j + 8] == 'A') || (arry[j + 6] == 'C' && arry[j + 7] == 'U' && arry[j + 8] == 'A'))){
                //printf("2's if! ");
                if (!((arry[j + 12] == 'U' && arry[j + 13] == 'A' && arry[j + 14] == 'A') || (arry[j + 12] == 'U' && arry[j + 13] == 'A' && arry[j + 14] == 'G')
                      || (arry[j + 12] == 'U' && arry[j + 13] == 'G' && arry[j + 14] == 'A'))){
                    //printf("3's if! ");
                    if (!((arry[j + 9] == 'U' && arry[j + 10] == 'A' && arry[j + 11] == 'A') || (arry[j + 9] == 'U' && arry[j + 10] == 'A' && arry[j + 11] == 'G')
                          || (arry[j + 9] == 'U' && arry[j + 10] == 'G' && arry[j + 11] == 'A') || (arry[j + 9] == 'U' && arry[j + 10] == 'C' && arry[j + 11] == 'A')
                          || (arry[j + 9] == 'U' && arry[j + 10] == 'U' && arry[j + 11] == 'A') || (arry[j + 9] == 'C' && arry[j + 10] == 'U' && arry[j + 11] == 'A'))){
                        //printf("4's if! ");
                        if (!(arry[j + 3] == 'U' && arry[j + 4] == 'A' && arry[j + 5] == 'C')){
                            //printf("5's if! ");
                            count = 1;// makes count 1
                            count2 = 0;
                            count3 = 0;
                            for (k = 0; k < 30; k++){
                                temp[k] = arry[j+k]; //store 30bp of sequence in temp
                            }
                            // see if restriction site string is in temp string
                            for (k = 0; k < 30; k++){
                                //printf("result here:\n");
                                //printf("%c", arry[j + k]);//print the 30bp target RNA part
                            }
                            printf("\n");//newline
                            strcat(temp2, temp);//add target RNA sequence to switch
                            if (temp2[29] == 'A') {
                                strcat(temp2, uaurbs);//add uau rbs sequence if end with A to switch
                            }
                            else{
                                strcat(temp2, auarbs);//add aua rbs sequence if not end with A to switch
                            }
                            reversestring(temp);
                            complement(temp);//reverse complement target RNA sequence
                            temp[16] = '\0';//only 15bp of reverse is used
                            temp[3] = 'A';
                            temp[4] = 'U';
                            temp[5] = 'G';//add aug start codon
                            strcat(temp2, temp);//add the reverse complement to switch
                            strcpy(temp3, temp2);//copy switch RNA to switch DNA temporary storage
                            reversestring(temp3);//reverse complement the RNA and (U to T) to make it DNA
                            complement(temp3);
                            makedna(temp3);
                            if ((strstr(temp3, ecor1) == NULL) || (strstr(temp3, not1) == NULL) || (strstr(temp3, xba1) == NULL)){//check if restriction site is inside the sequence in the middle part(temp2)
                                strcat(temp2, link);//add end cutting site to switch
                                strcpy(temp4, temp2);//add the middle part
                                memset(temp3,'\0',strlen(temp3));//clear temp3 contents
                                strcpy(temp3, temp4);//copy switch RNA to switch DNA temporary storage
                                strcpy(seq, temp4);
                                strcpy(dimer, arry);//Trigger&Switch
                                strcat(dimer, "&");
                                strcat(dimer, seq);//Trigger&Switch

                                makedna(temp3);
                                
                                
                                
                                char  *mfe_structure = vrna_alloc(sizeof(char) * (strlen(seq) + 1));
                                char  *mfe_structure_dimer = vrna_alloc(sizeof(char) * (strlen(dimer) + 1));
                                /* get a vrna_fold_compound with MFE and PF DP matrices and default model details */
                                vrna_fold_compound_t *vc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
                                vrna_fold_compound_t *vc_dimer = vrna_fold_compound(dimer, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
                                /* call MFE function */
                                double mfe = (double)vrna_mfe(vc, mfe_structure);
                                double mfe_dimer = (double)vrna_mfe_dimer(vc_dimer, mfe_structure_dimer);
                                
                                //start making rbsdimer
                                
                                strcpy(rbsdimer, rbs);
                                for (x=0; x<strlen(rbsdimer); x++) {
                                    rbsdimer[x] = mfe_structure_dimer[l+33+x];//copies rbs in dimer location to rbsdimer
                                }
                                
                                if (choice == '0' || choice == '2') {
                                    for (k = 0; k < strlen(rbsdimer); k++){//make all '.' if chosen not to check
                                        rbsdimer[k] = '.';
                                    }
                                }

                                //finish making rbsdimer
                                
                                //start making AUGdimer
                                
                                strcpy(AUGdimer, AUG);
                                for (x=0; x<3; x++) {
                                    AUGdimer[x] = mfe_structure_dimer[l+54+x];//copies rbs in dimer location to rbsdimer
                                }
                                    
                                if (choice == '0' || choice == '1') {
                                    for (k = 0; k < strlen(AUGdimer); k++){//make all '.' if chosen not to check
                                        AUGdimer[k] = '.';
                                    }
                                }
                                
                                //finish making AUGdimer
                                
                                printf("RBS%s dimerRBS%s AUG%s dimerAUG%s\n", rbstemp, rbsdimer, AUG, AUGdimer);
                                if ((strcmp(rbstemp, rbsdimer) == 0) && (strcmp(AUG, AUGdimer) == 0)) {//if rbs in dimer AND AUG in dimer is all '.'

                                    correctfolding = 1;
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                //start comparing mfe
                                if (newwindow == 1) {//clear data for new window comparison
                                    mfedifference = DBL_MIN;
                                    templocation = 0;
                                    for (k = 0; k < strlen(bestswitch); k++){//clear data
                                        bestswitch[k] = '\0';
                                    }
                                }
                                if ((mfe-mfe_dimer) > mfedifference) {
                                    mfedifference = (mfe-mfe_dimer);
                                    mfetemp = mfe;
                                    templocation = location;
                                    strcpy(bestswitch, seq);
                                }
                                
                                //end comparing mfe
                                
                                fprintf(textfile, "Group number: %d\nLocation: %d\nSequence: ", window, location);
                                printf("\nGroup number:%d \nDNA sequence: ", window);// print in txt file
                                for (k = 0; k < strlen(temp3); k++){//print final switch sequence
                                    printf("%c", temp3[k]);
                                }
                                for (k = 0; k < strlen(temp4); k++){//print final switch sequence in txt file
                                    fprintf(textfile, "%c", temp4[k]);
                                }
                                printf("\nRNA sequence: ");
                                for (k = 0; k < strlen(temp4); k++){//print final switch sequence
                                    printf("%c", temp4[k]);
                                }
                                fprintf(textfile, "\nDNA Sequence: ");// print in txt file
                                for (k = 0; k < strlen(temp3); k++){//print final switch sequence in txt file
                                    fprintf(textfile, "%c", temp3[k]);
                                }
                                fprintf(textfile, "\nMFE: %fkcal/mol\nMFE Structure: %s\nMFE Dimer: %fkcal/mol\nMFE Dimer Structure: %s\n", mfe, mfe_structure, mfe_dimer, mfe_structure_dimer);
                                fprintf(textfile, "\n\n");// print in txt file
                                printf("\nMFE: %fkcal/mol\nMFE Structure: %s\nMFE Dimer: %fkcal/mol\nMFE Dimer Structure: %s\n", mfe, mfe_structure, mfe_dimer, mfe_structure_dimer);
                                printf("\n\n");
                                    
                                    
                                    
                                    
                                }//end if(rbs and aug are all '.'s)

                                
                                
                                
                                
                                
                                
                                
                                
                                for (k = 0; k < strlen(temp); k++){//clear data
                                    temp[k] = '\0';
                                }
                                for (k = 0; k < strlen(temp2); k++){//clear data
                                    temp2[k] = '\0';
                                }
                                for (k = 0; k < strlen(temp3); k++){//clear data
                                    temp3[k] = '\0';
                                }
                                for (k = 0; k < strlen(temp4); k++){//clear data
                                    temp4[k] = '\0';
                                }
                                for (k = 0; k < strlen(seq); k++){//clear data
                                    seq[k] = '\0';
                                }
                                for (k = 0; k < strlen(rbsdimer); k++){//clear data
                                    rbsdimer[k] = '\0';
                                }
                                for (k = 0; k < strlen(AUGdimer); k++){//clear data
                                    AUGdimer[k] = '\0';
                                }
                                /* free mfe structure */
                                free(mfe_structure);
                                free(mfe_structure_dimer);
                                /* free memory occupied by vrna_fold_compound */
                                vrna_fold_compound_free(vc);
                                vrna_fold_compound_free(vc_dimer);

                            }

                            
                            
                            
                            
                            
                            
                            
                            //
                            //printf("finished");
                        }
                    }
                }
            }
        }
    }
                                return 0;
}



int main(void) {

    int window = 0;
    char input[len] = {'0'};
    char rbs[500] = {'0'};
    int repeat = 0;
    char copy[len] = {'0'};
    int i, j, l;
    int n = 0;
    int between;
    int choice;//RBS and AUG choice
    FILE *textfile;
    //create data.txt file then writes output in
    printf("Hello\n");
    textfile = fopen("data1.txt", "w");
    fprintf(textfile, "Hello!\n\n\n\n");
    
    do {
        printf("Enter RNA sequence: ");
        fgets(input, 1001, stdin);
        strtok(input, "\n");//removes '\n' from input
        makerna(input);
    }while (check_valid(input) != 0);
    do {
        printf("Enter RBS sequence(Type 'enter' for standard rbs{AUUAAAGAGGAGAAA}): ");
        fgets(rbs, 1001, stdin);
        if (rbs[0] == '\n') {
            strcpy(rbs, "AUUAAAGAGGAGAAA");
        }
        makerna(rbs);
        strtok(rbs, "\n");//removes '\n' from rbs
        } while (check_valid(rbs) != 0);
    do {
        printf("Enter distance between groups: ");
        scanf("%d", &between);
        printf("Eliminate for RBS and/or AUG pairings? ('0' for nocheck, '1' for RBS, '2' for AUG,'3' for BOTH): ");
        scanf("%d", &choice);
        //printf("\n get\n");
    } while (choice_valid(choice) != 0);
    
    l = strlen(input);
    strcpy(copy, input);
    
    //printf("123\n");

    look(choice, between, copy, rbs, l, window, n, textfile);

    //printf("%d\n",round);
    


    fclose(textfile);
    
    printf("-----------------------------------------\n");
    return 0;
}
