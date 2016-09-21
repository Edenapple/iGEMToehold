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

#include  "data_structures.h"
#include  "params.h"
#include  "utils.h"
#include  "eval.h"
#include  "fold.h"
#include  "part_func.h"
#include  "plot_structure.h"

#define len  100001



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

void complementdna(char *string) {// make complement
    int i;
    for (i=0; i<strlen(string); i++) {
        if (string[i] == 'A')
            string[i] = 'T';
        else if (string[i] == 'C')
            string[i] = 'G';
        else if (string[i] == 'T')
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

void Upcase(char *str) {//unused
    int i;
    for (i = 0; i < strlen(str); i++) {
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

int look(char filename[len], double mfearray[len], char promoterseq[100], int triggerlength, double inputtemperature, int maxdomainpair, int maxrbspair, int promoter, int choice, int printstructure, int plotrnachoice, int between, char res5[100], char res3[100], char arry[len], char rbs[500], int l, int window, int n, FILE *textfile, FILE *csvfile) {
    char RNAplot[1000] = {'\0'};//name of poscript plot file
    char switchnumberarray[100] = {'\0'};//switch number from int to char
    char rbslinker[1000] = {'\0'};//rbslinker
    char trigger[1000] = {'\0'};//trigger
    char trigger1[1000] = {'\0'};//trigger1
    char trigger30[10000] = {'\0'};//30bp trigger
    char trigger30temp[1000] = {'\0'};//30bp triggertemp
    char ecor1[] = "GAATTC"; //cutting site
    char not1[] = "GCGGCCGC"; //cutting site
    char xba1[] = "TCTAGA"; //cutting site
    char uaurbs[1000] = "UAU"; //rbs site, starts with UAU
    char auarbs[1000] = "AUA"; //rbs site, starts with AUA
    char addrestriction[len] = {'\0'};//switch with restriction
    char link[] = "AACCUGGCGGCAGCGCAAAAG";//linker sequence + end cutting site
    char temp[31] = {'\0'};//store the 30bp RNA thing([31] to store the null character)
    char temp2[1000] = {'\0'};//store toehold switch sequence
    char temp3[1000] = {'\0'};//store toehold switch DNA sequence
    char temp4[1000] = {'\0'};//temporary storage for additional things
    char temp5[1000] = {'\0'};//extra RBS thing
    char temp6[1000] = {'\0'};//reverse temp3 for restriction checking
    char seq[1000] = {'\0'};//temporary storage free energy process
    char dimer[10000] = {'\0'};//trigger switch dimer sequence
    char dimerswitch[10000] = {'\0'};//trigger switch dimer sequence
    char dimertrigger[10000] = {'\0'};//trigger switch dimer sequence
    char rbstemp[500] = {'\0'};;//turns all rbs things into '.'
    char rbsdimer[500] = {'\0'};;//rbs inside dimer
    char AUG[4] = "...";//chekc if AUG is "..."
    char AUGdimer[4] = "000";//AUG in dimer
    double mfearray2[len] = {'\0'};//store mfe
    char switchprimerforward[len] = {'\0'};
    char switchprimerreverse[len] = {'\0'};
    char triggerprimerforward[len] = {'\0'};
    char triggerprimerreverse[len] = {'\0'};
    int domainpair = 0;
    int domainunpair = 0;
    int rbspair = 0;
    int rbsunpair = 0;
    int rbslinkerpair = 0;
    int rbslinkerunpair = 0;
    int a = 0;//if T7 add ggg, then all calculations must +3
    int j;
    int k;
    int x;
    int y;
    int z=0;
    int count = 0;// if count = 0, window++, successful switch makes count 1, unsuccessful makes 2;
    int count2 = 0;// to prevent window gain after consecutive unsuccessful switches
    int count3 = 0;//space between each group
    int location = 0;// location
    int newwindow = 2;//new window or not
    int correctfolding = 0;//correct dimer folding rbs and aug
    double mfedifference = DBL_MIN;//difference between mfe for comparison
    double mfetemp = 0;//mfe of bestswitch
    char bestswitch[1000] = {'\0'};//best switch with lowest mfe in a group
    int templocation = 0;//location of best switch
    int switchnumber = 0;//switch number
    
    char triggertemp[1000] = {'\0'};
    char dimerstructuretemp[len] = {'\0'};
    double mfedimertemp = DBL_MAX;
    
    
    printf("\n\nstarted LOOKING\n\n\n");
    
    
    
    vrna_md_defaults_temperature(inputtemperature);

    
    
    
    a = 0;
    if (promoter == 1) {
        a = 3;//T7 or not
    }
    
    strcpy(rbstemp, rbs);
    for (j=0; j<strlen(rbstemp); j++) {
        rbstemp[j] = '.';//turns all rbstemp things into '.'
    }
    
    strcat(uaurbs, rbs);
    strcat(uaurbs, "AUA");
    strcat(auarbs, rbs);
    strcat(auarbs, "UAU");
    
    if (printstructure == 1) {
    fprintf(csvfile, "Switch Number,Group Number,Location,RNA Sequence,DNA Sequence,%dbp Trigger Sequence Location,%dbp Trigger Sequence,30bp Trigger Sequence,Switch Primer Forward,Switch Primer Reverse,Trigger Primer Forward,Trigger Primer Reverse,Overlapping Sequence for Switch Primer,Overlapping Sequence for Trigger Primer,Number of Toehold Domain Paired Bases,Number of Toehold Domain Unpaired Bases,Number of RBS Domain Paired Bases in %dbp Trigger-Switch Dimer,Number of RBS Domain Unpaired Bases in %dbp Trigger-Switch Dimer,MFE Switch,MFE  Switch Structure,MFE %dbp Trigger,MFE %dbp Trigger Structure,MFE Switch-Switch Dimer,MFE Switch-Switch Dimer Structure,MFE %dbp Trigger-Trigger Dimer,MFE %dbp Trigger-Trigger Dimer Structure,MFE %dbp Trigger-Switch Dimer,MFE %dbp Trigger-Switch Dimer Structure,MFE 30bp Trigger-Switch Dimer,MFE 30bp Trigger-Switch Dimer Structure,MFE RBS-Linker,MFE RBS-Linker Structure,Number of RBS-Linker Domain Paired Bases,Number of RBS-Linker Domain Unpaired Bases,MFE Difference,MFE Difference * MFE Switch\n", triggerlength,  triggerlength, triggerlength, triggerlength, triggerlength, triggerlength, triggerlength, triggerlength, triggerlength, triggerlength);
        }
    else{
        fprintf(csvfile, "Switch Number,Group Number,Location,RNA Sequence,DNA Sequence,%dbp Trigger Sequence Location,%dbp Trigger Sequence,30bp Trigger Sequence,Switch Primer Forward,Switch Primer Reverse,Trigger Primer Forward,Trigger Primer Reverse,Overlapping Sequence for Switch Primer,Overlapping Sequence for Trigger Primer,Number of Toehold Domain Paired Bases,Number of Toehold Domain Unpaired Bases,Number of RBS Domain Paired Bases in %dbp Trigger-Switch Dimer,Number of RBS Domain Unpaired Bases in %dbp Trigger-Switch Dimer,MFE Switch,MFE %dbp Trigger,MFE Switch-Switch Dimer,MFE %dbp Trigger-Trigger Dimer,MFE %dbp Trigger-Switch Dimer,MFE 30bp Trigger-Switch Dimer,MFE RBS-Linker,Number of RBS-Linker Domain Paired Bases,Number of RBS-Linker Domain Unpaired Bases,MFE Difference,MFE Difference * MFE Switch\n", triggerlength,  triggerlength, triggerlength, triggerlength, triggerlength, triggerlength, triggerlength);
    }
    
    if (l < 31){
        printf("Too short! ");
        return 0;
    }
    for (j = 0; j<(strlen(arry)-(triggerlength-1)); j++) {
        
        printf("\n\nstarted LOOKING(%d)\n\n\n", j);
        
        //reset strings in case illegal switch(other reset things is inside if statement)
        memset(temp,'\0',strlen(temp));
        memset(temp2,'\0',strlen(temp2));
        memset(temp3,'\0',strlen(temp3));
        memset(temp4,'\0',strlen(temp4));
        memset(temp5,'\0',strlen(temp5));
        memset(temp6,'\0',strlen(temp6));
        memset(seq,'\0',strlen(seq));
        memset(rbsdimer,'\0',strlen(rbsdimer));
        memset(AUGdimer,'\0',strlen(AUGdimer));
        memset(trigger,'\0',strlen(trigger));
        memset(trigger1,'\0',strlen(trigger1));
        memset(addrestriction,'\0',strlen(addrestriction));
        memset(rbslinker,'\0',strlen(rbslinker));
        
        if (count == 0){
            if (count2 == 0 && count3 == between) {
                newwindow = 1;//newwindow true
                window++;//new window
                count3 = 0;
                count2 = 1;
                if((window>0) && (correctfolding == 1)){
                    correctfolding = 0;
                //start printing best switch of the group
                fprintf(textfile, "\n\nSwitch with highest mfe difference in group number %d\nLocation: %d\nSequence: ", (window-1), templocation);
                    for (k = 0; k < strlen(bestswitch); k++){//print final switch sequence in txt file
                        fprintf(textfile, "%c", bestswitch[k]);
                    }
                    makedna(bestswitch);
                fprintf(textfile, "\nDNA Sequence: ");
                    for (k = 0; k < strlen(bestswitch); k++){//print final switch sequence in txt file
                        fprintf(textfile, "%c", bestswitch[k]);
                    }
                fprintf(textfile, "\nMFE: %fkcal/mol\nMFE difference: %fkcal/mol\n\n\n\n\n\n\n", mfetemp, mfedifference);
                //end printing best switch of the
                    mfedifference = DBL_MIN;
                }
            }
            else{
                newwindow = 0;//newwindow false
            }
            count = 2;
        }
        
        
        
        if ((j)>((int)(strlen(arry))-triggerlength)) {
            break;
        }
        
        
        
        
        printf("\nIt works1\n");
        
        
        count = 0;//reset count
        count3++;
        location=j;//location increase
        if ((arry[j + 12] != 'A' && arry[j + 12] != 'U' && arry[j + 13] != 'A' && arry[j + 13] != 'U' && arry[j + 14] != 'A' && arry[j + 14] != 'U') || (arry[j + 12] != 'A' && arry[j + 12] != 'U' && arry[j + 13] != 'A' && arry[j + 13] != 'U') || (arry[j + 13] != 'A' && arry[j + 13] != 'U' && arry[j + 14] != 'A' && arry[j + 14] != 'U') || (arry[j + 12] != 'A' && arry[j + 12] != 'U' && arry[j + 14] != 'A' && arry[j + 14] != 'U')) {//Check for (12, 13, 14) has two consecutive CG base pair followed by an AU base pair (i.e. sequence is ACC, ACG, UCC, UCG, AGC, UGC, AGG,UGG). If it has=> valid (if final output is zero, count the number of G/C for these three bp, if it >=2, =>valid)
            printf("\n1's if! ");
            if (!((arry[j + 6] == 'U' && arry[j + 7] == 'A' && arry[j + 8] == 'A') || (arry[j + 6] == 'U' && arry[j + 7] == 'A' && arry[j + 8] == 'G')
                  || (arry[j + 6] == 'U' && arry[j + 7] == 'G' && arry[j + 8] == 'A') || (arry[j + 6] == 'U' && arry[j + 7] == 'C' && arry[j + 8] == 'A')
                  || (arry[j + 6] == 'U' && arry[j + 7] == 'U' && arry[j + 8] == 'A') || (arry[j + 6] == 'C' && arry[j + 7] == 'U' && arry[j + 8] == 'A'))){
                printf("\n2's if! ");
                if (!((arry[j + 12] == 'U' && arry[j + 13] == 'A' && arry[j + 14] == 'A') || (arry[j + 12] == 'U' && arry[j + 13] == 'A' && arry[j + 14] == 'G')
                      || (arry[j + 12] == 'U' && arry[j + 13] == 'G' && arry[j + 14] == 'A'))){
                    printf("\n3's if! ");
                    if (!((arry[j + 9] == 'U' && arry[j + 10] == 'A' && arry[j + 11] == 'A') || (arry[j + 9] == 'U' && arry[j + 10] == 'A' && arry[j + 11] == 'G')
                          || (arry[j + 9] == 'U' && arry[j + 10] == 'G' && arry[j + 11] == 'A') || (arry[j + 9] == 'U' && arry[j + 10] == 'C' && arry[j + 11] == 'A')
                          || (arry[j + 9] == 'U' && arry[j + 10] == 'U' && arry[j + 11] == 'A') || (arry[j + 9] == 'C' && arry[j + 10] == 'U' && arry[j + 11] == 'A'))){
                        printf("\n4's if! ");
                        if (!(arry[j + 3] == 'U' && arry[j + 4] == 'A' && arry[j + 5] == 'C')){
                            printf("\n5's if! ");
                            
                            
                            printf("\nIt works2\n");
                            for (k = 0; k < 30; k++){
                                temp[k] = arry[j+k]; //store 30bp of sequence in temp
                            }
                            // see if restriction site string is in temp string
                            
                            if (((strstr(temp, "AAAAA") == NULL) && (strstr(temp, "CCCCC") == NULL) && (strstr(temp, "UUUUU") == NULL) && (strstr(temp, "GGGGG") == NULL))){
                                
                                count = 1;// makes count 1
                                count2 = 0;
                                count3 = 0;
                                
                            strcpy(trigger30, temp);
                                strcpy(trigger30temp, temp);
                            
                                
                                reversestring(temp);//new
                                complement(temp);//new
                                
                            printf("\n");//newline
                                if (promoter == 1) {
                                    strcpy(temp2, "GGG");//GGG for T7
                                }
                                
                                
                            strcat(temp2, temp);//add target RNA sequence to switch
                                
                                
                            if (temp2[29+a] == 'A') {
                                strcat(temp2, uaurbs);//add uau rbs sequence if end with A to switch
                                strcpy(rbslinker, uaurbs);
                            }
                            else{
                                strcat(temp2, auarbs);//add aua rbs sequence if not end with A to switch
                                strcpy(rbslinker, auarbs);
                                strcat(temp5, uaurbs);
                            }
                            reversestring(temp);
                            complement(temp);//reverse complement target RNA sequence
                            temp[15] = '\0';//only 15bp of reverse is used
                            temp[3] = 'A';
                            temp[4] = 'U';
                            temp[5] = 'G';//add aug start codon
                            
                                
                            strcat(temp2, temp);//add the reverse complement to switch
                            strcat(rbslinker, temp);
                            strcat(temp5, temp);
                            strcpy(temp3, temp2);//copy switch RNA to switch DNA temporary storage
                            reversestring(temp3);//reverse complement the RNA and (U to T) to make it DNA
                            complement(temp3);
                            makedna(temp3);
                            reversestring(temp5);//same for temp5
                            complement(temp5);
                            makedna(temp5);
                            strcpy(temp6, temp3);
                            reversestring(temp6);
                            if ((strstr(temp3, ecor1) == NULL) && (strstr(temp3, not1) == NULL) && (strstr(temp3, xba1) == NULL) &&  (strstr(temp3, "CTCGAG") == NULL) &&  (strstr(temp3, "ACTAGT") == NULL) &&  (strstr(temp3, "CTGCAG") == NULL) && (strstr(temp6, ecor1) == NULL) && (strstr(temp6, not1) == NULL) && (strstr(temp6, xba1) == NULL) &&  (strstr(temp6, "CTCGAG") == NULL) &&  (strstr(temp6, "ACTAGT") == NULL) &&  (strstr(temp6, "CTGCAG") == NULL)){//check if restriction site is inside the sequence in the middle part(temp2)
                                
                                strcat(temp2, link);//add link to switch
                                strcat(rbslinker, link);
                                
                                strcat(addrestriction, temp2);
                                strcat(addrestriction, res3);//3 end restriction
                                
                                strcpy(temp4, addrestriction);//copy addrestricito to temp4
                                memset(temp3,'\0',strlen(temp3));//clear temp3 contents
                                strcpy(temp3, temp4);//copy switch RNA(temp4) to switch DNA temporary storage(temp3)
                                strcpy(seq, temp4);
                                printf("Location: %d\n\n", location);
                                /*start making trigger
                                printf("Input-100:%d\n", (int)(strlen(arry))-100);
                                if ((location)<((int)(strlen(arry))-100)) {
                                    for (k = 0; k < 100; k++){//clear data
                                        trigger1[k] = arry[k+(location)];
                                    }
                                }
                                else{
                                for (k = 0; k < 100; k++){//clear data
                                    trigger1[k] = arry[k+((int)(strlen(arry))-100)];
                                }
                                }

                                
                                //strcat(trigger, trigger1);
                                //strcat(trigger, res3);//3 end restriction
                                //finish making trigger
                                */
                                
                                
                                //strcpy(dimer, trigger);//Trigger&Switch
                                //strcat(dimer, "&");
                                //strcat(dimer, seq);//Trigger&Switch
                                
                                
                                strcat(trigger30, "&");
                                strcat(trigger30, seq);//Trigger&Switch
                                
                                strcpy(dimerswitch, seq);//Switch&Switch
                                strcat(dimerswitch, "&");
                                strcat(dimerswitch, seq);//Switch&Switch
                                


                                makedna(temp3);
                                
                                printf("\nRBSLINKER: %s\n",rbslinker);
                                
                                char  *mfe_structure = vrna_alloc(sizeof(char) * (strlen(seq) + 1));
                                char  *mfe_structure_dimerswitch = vrna_alloc(sizeof(char) * (strlen(dimerswitch) + 1));
                                char  *mfe_structure_trigger30 = vrna_alloc(sizeof(char) * (strlen(trigger30) + 1));
                                char  *mfe_structure_rbslinker = vrna_alloc(sizeof(char) * (strlen(rbslinker) + 1));
                                /* get a vrna_fold_compound with MFE and PF DP matrices and default model details */
                                vrna_fold_compound_t *vc = vrna_fold_compound(seq, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
                                vrna_fold_compound_t *vc_dimerswitch = vrna_fold_compound(dimerswitch, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
                                vrna_fold_compound_t *vc_trigger30 = vrna_fold_compound(trigger30, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
                                vrna_fold_compound_t *vc_rbslinker = vrna_fold_compound(rbslinker, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
                                /* call MFE function */
                                double mfe = (double)vrna_mfe(vc, mfe_structure);
                                double mfe_dimerswitch = (double)vrna_mfe_dimer(vc_dimerswitch, mfe_structure_dimerswitch);
                                double mfe_trigger30 = (double)vrna_mfe_dimer(vc_trigger30, mfe_structure_trigger30);
                                double mfe_rbslinker = (double)vrna_mfe(vc_rbslinker, mfe_structure_rbslinker);
                                
                                
                                
                                
                                
                                //start calulating domain pair and domain unpair
                                
                                domainpair = 0;
                                domainunpair = 0;
                                
                                for (k=0; k<15+a; k++) {
                                    if (mfe_structure[k] == '.') {
                                        domainunpair++;
                                    }
                                    else{
                                        domainpair++;
                                    }
                                }
                                
                                //finish calulating domain pair and domain unpair
                                
                                if(domainpair<=maxdomainpair){
                                
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                mfedimertemp = DBL_MAX;
                                int rbspairtemp = INT_MAX;//for comparison
                                for (y=0; y<1/*45*/; y++) {
                                    //if ((location+y)<((int)(strlen(arry))-74)) {
                                        
                                        memset(trigger1,'\0',strlen(trigger1));
                                        memset(trigger,'\0',strlen(trigger));
                                        memset(dimer,'\0',strlen(dimer));
                                        
                                        for (k = 0; k < triggerlength; k++){//clear data
                                            trigger1[k] = arry[k+(location+y)];;
                                        }
                                    //}
                                    //else{
                                        //break;
                                    //}
                                    
                                    if (promoter == 1) {
                                        strcpy(trigger, "GGG");//T7
                                    }
                                    
                                    strcat(trigger, trigger1);
                                    //strcat(trigger, res3);//3 end restriction
                                    //finish making trigger
                                    
                                    
                                    strcpy(dimer, trigger);//Trigger&Switch
                                    strcat(dimer, "&");
                                    strcat(dimer, seq);//Trigger&Switch
                                    
                                    char  *mfe_structure_dimertemp = vrna_alloc(sizeof(char) * (strlen(dimer) + 1));
                                    vrna_fold_compound_t *vc_dimertemp = vrna_fold_compound(dimer, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
                                    double mfe_dimertemp = (double)vrna_mfe_dimer(vc_dimertemp, mfe_structure_dimertemp);
                                    
                                    

                                    
                                    
                                    rbspair = 0;//rbs base pairs in 100bp trigger switch duplex
                                    rbsunpair = 0;
                                    
                                    
                                    for (x = 0; x<strlen(rbs); x++) {
                                        if (mfe_structure_dimertemp[triggerlength+a+33+x] == '.') {
                                            rbsunpair++;
                                        }
                                        else{
                                            rbspair++;
                                        }
                                    }
                                    
                                    
                                    //if (1==1/*rbspair<=maxrbspair*/) {
                                        
                                        if ((rbspairtemp)>(rbspair) || rbspairtemp == INT_MAX) {
                                            rbspairtemp = rbspair;
                                            //if ((mfe-mfe_dimertemp)>(mfe-mfedimertemp) || mfedimertemp == DBL_MAX) {
                                            z=y;
                                            memset(triggertemp,'\0',strlen(triggertemp));
                                            strcpy(triggertemp, trigger);
                                            memset(dimerstructuretemp,'\0',strlen(dimerstructuretemp));
                                            strcpy(dimerstructuretemp, mfe_structure_dimertemp);
                                            mfedimertemp = mfe_dimertemp;
                                            //printf("\ntriggertemp: %s\n", triggertemp);
                                            //}
                                        }
                                    //}
                                    
                                }
                                    
                                    if ((j)>((int)(strlen(arry))-triggerlength)) {
                                        break;
                                    }
                                
                                
                                memset(trigger,'\0',strlen(trigger));
                                strcpy(trigger, triggertemp);
                                printf("trigger: %s\n\n", trigger);
                                
                                    char  *mfe_structure_trigger = vrna_alloc(sizeof(char) * (strlen(trigger) + 1));
                                    vrna_fold_compound_t *vc_trigger = vrna_fold_compound(trigger, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
                                    
                                    
                                    
                                    strcpy(dimertrigger, trigger);//Trigger&Trigger
                                    strcat(dimertrigger, "&");
                                    strcat(dimertrigger, trigger);//Trigger&Trigger
                                
                                
                                
                               
                                    

                            char  *mfe_structure_dimertrigger = vrna_alloc(sizeof(char) * (strlen(dimertrigger) + 1));
                            vrna_fold_compound_t *vc_dimertrigger = vrna_fold_compound(dimertrigger, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
                            double mfe_dimertrigger = (double)vrna_mfe_dimer(vc_dimertrigger, mfe_structure_dimertrigger);
                            double mfe_trigger = (double)vrna_mfe(vc_trigger, mfe_structure_trigger);
                                    
                                    
                                    
                                
                                
                                
                                
                                //start making rbsdimer
                                
                                memset(rbsdimer,'\0',strlen(rbsdimer));
                                    
                                    
                                strcpy(rbsdimer, rbs);
                                for (x=0; x<strlen(rbsdimer); x++) {
                                    rbsdimer[x] = dimerstructuretemp[triggerlength+a+33+x];//copies rbs in dimer location to rbsdimer
                                }
                                
                                /*if (choice == 0 || choice == 2) {
                                    for (k = 0; k < 15; k++){//make all '.' if chosen not to check
                                        rbsdimer[k] = '.';
                                    }
                                }*/

                                //finish making rbsdimer
                                
                                //start making AUGdimer
                                
                                /*strcpy(AUGdimer, AUG);
                                for (x=0; x<3; x++) {
                                    AUGdimer[x] = dimerstructuretemp[100+a+54+x];//copies rbs in dimer location to rbsdimer
                                }
                                    
                                if (choice == 0 || choice == 1) {
                                    for (k = 0; k < 3; k++){//make all '.' if chosen not to check
                                        AUGdimer[k] = '.';
                                    }
                                }*/
                                
                                //finish making AUGdimer
                                
                                    rbspair = 0;//rbs base pairs in 100bp trigger switch duplex
                                    rbsunpair = 0;
                                    
                                    
                                    for (x = 0; x<strlen(rbs); x++) {
                                        if (dimerstructuretemp[triggerlength+a+33+x] == '.') {
                                            rbsunpair++;
                                        }
                                        else{
                                            rbspair++;
                                        }
                                    }
                                    
                                    
                                    
                                    rbslinkerpair = 0;//rbslinker base pairs in 100bp trigger switch duplex
                                    rbslinkerunpair = 0;

                                    
                                    for (x = 0; x<strlen(rbslinker); x++) {
                                        if (dimerstructuretemp[x+a+triggerlength+33] == '.') {
                                            rbslinkerunpair++;
                                        }
                                        else{
                                            rbslinkerpair++;
                                        }
                                    }
                                    
                                printf("Choice %d RBS%s dimerRBS%s AUG%s dimerAUG%s\n", choice, rbstemp, rbsdimer, AUG, AUGdimer);
                                    
                                    
                                    
                                //if ((strcmp(rbstemp, rbsdimer) == 0) && (strcmp(AUG, AUGdimer) == 0)) {//if rbs in dimer AND AUG in dimer is all '.'}
                                    if (rbspair<=maxrbspair) {
                                    

                                    correctfolding = 1;
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                //start comparing mfe
                                if (newwindow == 1) {//clear data for new window comparison
                                    mfedifference = DBL_MIN;
                                    templocation = 0;
                                    for (k = 0; k < strlen(bestswitch); k++){//clear data
                                        bestswitch[k] = '\0';
                                    }
                                }
                                if ((mfe-mfedimertemp) > mfedifference) {
                                    mfedifference = (mfe-mfedimertemp);
                                    mfetemp = mfe;
                                    templocation = location;
                                    strcpy(bestswitch, seq);
                                }
                                
                                //end comparing mfe
                                
                                    mfearray[switchnumber+1] = mfe-mfedimertemp;//store mfe difference
                                    
                                    switchnumber++;
                                    
                                fprintf(textfile, "Switch Number: %d\nGroup number: %d\nLocation: %d\nSequence: ", switchnumber, window, location);
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
                                fprintf(textfile, "\n150bp Trigger Sequence: ");// print in txt file
                                for (k = 0; k < strlen(trigger); k++){//print final switch sequence in txt file
                                    fprintf(textfile, "%c", trigger[k]);
                                }
                                fprintf(textfile, "\n30bp Trigger Sequence: ");// print in txt file
                                for (k = 0; k < 30; k++){//print final switch sequence in txt file
                                    fprintf(textfile, "%c", trigger30[k]);
                                }
                                fprintf(textfile, "\nNumber of Toehold Domain Paired Bases: %d\nNumber of Toehold Domain Unpaired Bases: %d\nMFE: %fkcal/mol\nMFE Structure: %s\nMFE Trigger: %fkcal/mol\nMFE Trigger Structure: %s\nMFE Switch Dimer: %fkcal/mol\nMFE Switch Dimer Structure: %s\nMFE Trigger Self Dimer: %fkcal/mol\nMFE Trigger Self Dimer Structure: %s\nMFE 150bp Trigger Dimer: %fkcal/mol\nMFE 150bp Trigger Dimer Structure: %s\nMFE 30bp Trigger Dimer: %fkcal/mol\nMFE 30bp Trigger Dimer Structure: %s\nMFE RBS Linker: %f\nMFE RBS Linker Structure: %s\n", domainpair, domainunpair, mfe, mfe_structure, mfe_trigger, mfe_structure_trigger, mfe_dimerswitch, mfe_structure_dimerswitch, mfe_dimertrigger, mfe_structure_dimertrigger, mfedimertemp, dimerstructuretemp, mfe_trigger30, mfe_structure_trigger30, mfe_rbslinker, mfe_structure_rbslinker);
                                fprintf(textfile, "\n\n");// print in txt file
                                printf("\nMFE: %fkcal/mol\nMFE Structure: %s\nMFE Dimer: %fkcal/mol\nMFE Dimer Structure: %s\n", mfe, mfe_structure, mfedimertemp, dimerstructuretemp);
                                printf("\n\n");
                                 
                                
                                    
                                    
                                    memset(temp6,'\0',strlen(temp6));
                                    strcpy(temp6, res5);
                                    strcat(temp6, promoterseq);
                                    strcat(temp6, seq);
                                        
                                    memset(temp3,'\0',strlen(temp3));
                                    strcpy(temp3, temp6);
                                    makedna(temp3);
                                    
                                    
                                    memset(temp6,'\0',strlen(temp6));
                                    strcpy(temp6, res5);
                                    strcat(temp6, promoterseq);
                                    strcat(temp6, trigger);
                                    strcat(temp6, res3);
                                    memset(trigger,'\0',strlen(trigger));
                                    strcpy(trigger, temp6);
                                    
                                        
                                        
                                    //start making primers
                                        int CGcountfoward = 0;
                                        int CGcountreverse = 0;
                                        int overlap = 30;
                                        char switchoverlap[len] = {'\0'};
                                        char triggeroverlap[len] = {'\0'};
                                        memset(switchprimerforward,'\0',strlen(switchprimerforward));
                                        memset(switchprimerreverse,'\0',strlen(switchprimerreverse));
                                        memset(triggerprimerforward,'\0',strlen(triggerprimerforward));
                                        memset(triggerprimerreverse,'\0',strlen(triggerprimerreverse));
                                        
                                        /*
                                        for (x=0; x<6; x++) {
                                            overlap = 2*(15-x);
                                            memset(switchprimerforward,'\0',strlen(switchprimerforward));
                                            memset(switchprimerreverse,'\0',strlen(switchprimerreverse));
                                            for (k=0; k<(((strlen(temp3))/2)+15-x); k++) {
                                                switchprimerforward[k] = temp3[k];
                                                switchprimerreverse[k] = temp3[strlen(temp3)-1-k];
                                            }
                                            CGcountfoward = 0;
                                            CGcountreverse = 0;
                                            for (y=0; y<strlen(switchprimerforward); y++) {
                                                if (switchprimerforward[y] == 'C' || switchprimerforward[y] == 'G') {
                                                    CGcountfoward++;
                                                }
                                            }
                                            for (y=0; y<strlen(switchprimerreverse); y++) {
                                                if (switchprimerreverse[y] == 'C' || switchprimerreverse[y] == 'G') {
                                                    CGcountreverse++;
                                                }
                                            }
                                            if ((double)CGcountfoward/(double)strlen(switchprimerforward) > 0.4 && (double)CGcountfoward/(double)strlen(switchprimerforward) < 0.6 && (double)CGcountreverse/(double)strlen(switchprimerreverse) > 0.4 && (double)CGcountreverse/(double)strlen(switchprimerreverse) < 0.6) {
                                                break;
                                            }
                                        }
                                        
                                        for (y=0; y<overlap; y++) {
                                            switchoverlap[y] = temp3[y+strlen(switchprimerforward)-(overlap/2)];
                                        }
                                        
                                        */
                                        
                                        memset(switchprimerforward,'\0',strlen(switchprimerforward));
                                        memset(switchprimerreverse,'\0',strlen(switchprimerreverse));
                                        for (k=0; k<(((strlen(temp3))/2)+15); k++) {
                                            switchprimerforward[k] = temp3[k];
                                            switchprimerreverse[k] = temp3[strlen(temp3)-1-k];
                                        }
                                        
                                        memset(switchoverlap,'\0',strlen(switchoverlap));
                                        for (y=0; y<overlap; y++) {
                                            switchoverlap[y] = temp3[y+strlen(switchprimerforward)-(overlap)];
                                        }
                                        
                                        printf("\nSwitchoverlap: %s\n", switchoverlap);
                                        
                                        CGcountfoward = 0;
                                        CGcountreverse = 0;
                                        for (y=0; y<strlen(switchprimerforward); y++) {
                                            if (switchprimerforward[y] == 'C' || switchprimerforward[y] == 'G') {
                                                CGcountfoward++;
                                            }
                                        }
                                        for (y=0; y<strlen(switchprimerreverse); y++) {
                                            if (switchprimerreverse[y] == 'C' || switchprimerreverse[y] == 'G') {
                                                CGcountreverse++;
                                            }
                                        }

                                        
                                        x = strlen(switchprimerforward);
                                            for (k=0; ((double)CGcountfoward/(double)strlen(switchprimerforward) > 0.4 && (double)CGcountfoward/(double)strlen(switchprimerforward) < 0.6 && (double)CGcountreverse/(double)strlen(switchprimerreverse) > 0.4 && (double)CGcountreverse/(double)strlen(switchprimerreverse) < 0.6) && ((switchprimerforward[strlen(switchprimerforward)-1] == 'G') || (switchprimerforward[strlen(switchprimerforward)-1] == 'C')) && ((switchprimerreverse[strlen(switchprimerreverse)-1] == 'G') || (switchprimerreverse[strlen(switchprimerreverse)-1] == 'C')); k++) {
                                                memset(switchprimerforward,'\0',strlen(switchprimerforward));
                                                //printf("\nSwitchoverlap: %s\n", switchoverlap);
                                                for (y=0; y<(x-k); y++) {
                                                    switchprimerforward[y] = temp3[y];
                                                }
                                                memset(switchprimerreverse,'\0',strlen(switchprimerreverse));
                                                for (y=0; y<(x+k); y++) {
                                                    switchprimerreverse[y] = temp3[strlen(temp3)-1-y];
                                                }
                                                
                                                memset(switchoverlap,'\0',strlen(switchoverlap));
                                                for (y=0; y<overlap; y++) {
                                                    switchoverlap[y] = temp3[y+strlen(switchprimerforward)-(overlap)];
                                                }
                                                
                                                CGcountfoward = 0;
                                                CGcountreverse = 0;
                                                for (y=0; y<strlen(switchprimerforward); y++) {
                                                    if (switchprimerforward[y] == 'C' || switchprimerforward[y] == 'G') {
                                                        CGcountfoward++;
                                                    }
                                                }
                                                for (y=0; y<strlen(switchprimerreverse); y++) {
                                                    if (switchprimerreverse[y] == 'C' || switchprimerreverse[y] == 'G') {
                                                        CGcountreverse++;
                                                    }
                                                }
                                                
                                                
                                                makedna(switchprimerforward);
                                                makedna(switchprimerreverse);
                                                //reversestring(switchprimerreverse);
                                                complementdna(switchprimerreverse);
                                            }
                                        
                                    
                                        
                                        
                                        
                                        /*
                                        for (x=0; x<6; x++) {
                                            overlap = 2*(15-x);
                                            memset(triggerprimerforward,'\0',strlen(triggerprimerforward));
                                            memset(triggerprimerreverse,'\0',strlen(triggerprimerreverse));
                                            for (k=0; k<(((strlen(trigger))/2)+15-x); k++) {
                                                triggerprimerforward[k] = trigger[k];
                                                triggerprimerreverse[k] = trigger[strlen(trigger)-1-k];
                                            }
                                            CGcountfoward = 0;
                                            CGcountreverse = 0;
                                            for (y=0; y<strlen(triggerprimerforward); y++) {
                                                if (triggerprimerforward[y] == 'C' || triggerprimerforward[y] == 'G') {
                                                    CGcountfoward++;
                                                }
                                            }
                                            for (y=0; y<strlen(triggerprimerreverse); y++) {
                                                if (triggerprimerreverse[y] == 'C' || triggerprimerreverse[y] == 'G') {
                                                    CGcountreverse++;
                                                }
                                            }
                                            if ((double)CGcountfoward/(double)strlen(triggerprimerforward) > 0.4 && (double)CGcountfoward/(double)strlen(triggerprimerforward) < 0.6 && (double)CGcountreverse/(double)strlen(triggerprimerreverse) > 0.4 && (double)CGcountreverse/(double)strlen(triggerprimerreverse) < 0.6) {
                                                break;
                                            }
                                        }
                                        
                                        for (y=0; y<overlap; y++) {
                                            triggeroverlap[y] = trigger[y+strlen(triggerprimerforward)-(overlap/2)];
                                        }
                                        */
                                        
                                        memset(triggerprimerforward,'\0',strlen(triggerprimerforward));
                                        memset(triggerprimerreverse,'\0',strlen(triggerprimerreverse));
                                        for (k=0; k<(((strlen(trigger))/2)+15); k++) {
                                            triggerprimerforward[k] = trigger[k];
                                            triggerprimerreverse[k] = trigger[strlen(trigger)-1-k];
                                        }
                                        
                                        CGcountfoward = 0;
                                        CGcountreverse = 0;
                                        for (y=0; y<strlen(triggerprimerforward); y++) {
                                            if (triggerprimerforward[y] == 'C' || triggerprimerforward[y] == 'G') {
                                                CGcountfoward++;
                                            }
                                        }
                                        for (y=0; y<strlen(triggerprimerreverse); y++) {
                                            if (triggerprimerreverse[y] == 'C' || triggerprimerreverse[y] == 'G') {
                                                CGcountreverse++;
                                            }
                                        }
                                        
                                        memset(triggeroverlap,'\0',strlen(triggeroverlap));
                                        for (y=0; y<overlap; y++) {
                                            triggeroverlap[y] = trigger[y+strlen(triggerprimerforward)-(overlap)];
                                        }
                                        
                                        printf("\nTriggeroverlap: %s\n", triggeroverlap);
                                        
                                        x = strlen(triggerprimerforward);
                                        for (k=0; ((double)CGcountfoward/(double)strlen(triggerprimerforward) > 0.4 && (double)CGcountfoward/(double)strlen(triggerprimerforward) < 0.6 && (double)CGcountreverse/(double)strlen(triggerprimerreverse) > 0.4 && (double)CGcountreverse/(double)strlen(triggerprimerreverse) < 0.6) && ((triggerprimerforward[strlen(triggerprimerforward)-1] == 'G') || (triggerprimerforward[strlen(triggerprimerforward)-1] == 'C')) && ((triggerprimerforward[strlen(triggerprimerreverse)-1] == 'G') || (triggerprimerforward[strlen(triggerprimerreverse)-1] == 'C')); k++) {
                                            memset(triggerprimerforward,'\0',strlen(triggerprimerforward));
                                            memset(triggeroverlap,'\0',strlen(triggeroverlap));
                                            //printf("\nTriggeroverlap: %s\n", triggeroverlap);
                                            for (y=0; y<(x-k); y++) {
                                                triggerprimerforward[y] = trigger[y];
                                            }
                                            memset(triggerprimerreverse,'\0',strlen(triggerprimerreverse));
                                            for (y=0; y<(x+k); y++) {
                                                triggerprimerreverse[y] = trigger[strlen(trigger)-1-y];
                                            }
                                            
                                            for (y=0; y<overlap; y++) {
                                                triggeroverlap[y] = trigger[y+strlen(triggerprimerforward)-(overlap)];
                                            }
                                            
                                            CGcountfoward = 0;
                                            CGcountreverse = 0;
                                            for (y=0; y<strlen(triggerprimerforward); y++) {
                                                if (triggerprimerforward[y] == 'C' || triggerprimerforward[y] == 'G') {
                                                    CGcountfoward++;
                                                }
                                            }
                                            for (y=0; y<strlen(triggerprimerreverse); y++) {
                                                if (triggerprimerreverse[y] == 'C' || triggerprimerreverse[y] == 'G') {
                                                    CGcountreverse++;
                                                }
                                            }
                                            
                                            
                                            
                                            makedna(triggerprimerforward);
                                            makedna(triggerprimerreverse);
                                            //reversestring(triggerprimerreverse);
                                            complementdna(triggerprimerreverse);
                                        }
                                        
                                        
                                        
                                        
                                    //finish making primers
                                        
                                        
                                //start csv
                                    printf("\nSTART CSV\n");
                                        
                                        
                                    mfearray2[switchnumber] = mfe;//store mfe difference
                                        printf("\nmfearray2: %f\n", mfearray2[switchnumber]);
                                        if (printstructure == 1) {
                                            fprintf(csvfile, "%d,%d,%d,%s,%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,%f,%s,%f,%s,%f,%s,%f,%s,%f,%s,%f,%s,%f,%s,%d,%d,%f,%f\n", switchnumber, window, location, seq, temp3, z, trigger, trigger30temp, switchprimerforward, switchprimerreverse, triggerprimerforward, triggerprimerreverse, switchoverlap, triggeroverlap, domainpair, domainunpair, rbspair, rbsunpair, mfe, mfe_structure, mfe_trigger, mfe_structure_trigger, mfe_dimerswitch, mfe_structure_dimerswitch, mfe_dimertrigger, mfe_structure_dimertrigger, mfedimertemp, dimerstructuretemp, mfe_trigger30, mfe_structure_trigger30, mfe_rbslinker, mfe_structure_rbslinker, rbslinkerpair, rbslinkerunpair, (mfe+mfe_trigger-mfedimertemp), (mfe)*(mfe+mfe_trigger-mfedimertemp));
                                        }
                                        else{
                                            fprintf(csvfile, "%d,%d,%d,%s,%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f\n", switchnumber, window, location, seq, temp3, z, trigger, trigger30temp, switchprimerforward, switchprimerreverse, triggerprimerforward, triggerprimerreverse, switchoverlap, triggeroverlap, domainpair, domainunpair, rbspair, rbsunpair, mfe, mfe_trigger, mfe_dimerswitch, mfe_dimertrigger, mfedimertemp, mfe_trigger30, mfe_rbslinker, rbslinkerpair, rbslinkerunpair, (mfe+mfe_trigger-mfedimertemp), (mfe)*(mfe+mfe_trigger-mfedimertemp));
                                        }
                                    
                                    
                                    printf("\nEND CSV\n");
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                               //end csv
                                        if(plotrnachoice == 1){
                                    //start plot
                                    strcpy(RNAplot, "RNAplot ");
                                    sprintf(switchnumberarray, "%d", switchnumber);
                                    strcat(RNAplot, filename);
                                    strcat(RNAplot, " ");
                                    strcat(RNAplot, switchnumberarray);
                                    PS_rna_plot	(temp4, mfe_structure, RNAplot);
                                    
                                    memset(RNAplot,'\0',strlen(RNAplot));
                                    memset(switchnumberarray,'\0',strlen(switchnumberarray));
    
                                    //end plot
                                        }
                                    
                                }//end if(rbs and aug are all '.'s)

                                
                                
                                
                                
                                
                                
                                
                                //reset all strings
                                memset(temp,'\0',strlen(temp));
                                memset(temp2,'\0',strlen(temp2));
                                memset(temp3,'\0',strlen(temp3));
                                memset(temp4,'\0',strlen(temp4));
                                memset(temp5,'\0',strlen(temp5));
                                memset(temp6,'\0',strlen(temp6));
                                memset(seq,'\0',strlen(seq));
                                memset(dimer,'\0',strlen(dimer));
                                memset(dimerswitch,'\0',strlen(dimerswitch));
                                memset(dimertrigger,'\0',strlen(dimertrigger));
                                memset(rbsdimer,'\0',strlen(rbsdimer));
                                memset(AUGdimer,'\0',strlen(AUGdimer));
                                memset(trigger,'\0',strlen(trigger));
                                memset(trigger1,'\0',strlen(trigger1));
                                memset(trigger30,'\0',strlen(trigger1));
                                memset(addrestriction,'\0',strlen(addrestriction));
                                memset(rbslinker,'\0',strlen(rbslinker));
                                memset(dimerstructuretemp,'\0',strlen(dimerstructuretemp));

                                /* free mfe structure */
                                free(mfe_structure);
                                free(mfe_structure_trigger);
                                free(mfe_structure_trigger30);
                                free(mfe_structure_rbslinker);
                                free(mfe_structure_dimerswitch);
                                free(mfe_structure_dimertrigger);
                                /* free memory occupied by vrna_fold_compound */
                                vrna_fold_compound_free(vc);
                                vrna_fold_compound_free(vc_trigger);
                                vrna_fold_compound_free(vc_trigger30);
                                vrna_fold_compound_free(vc_rbslinker);
                                vrna_fold_compound_free(vc_dimerswitch);
                                vrna_fold_compound_free(vc_dimertrigger);

                            }

                            }
                            
                            }
                            
                            
                            
                            
                            //
                            //printf("finished");
                        }
                    }
                }
            }
        }
    }

    //
    window++;//new window
        correctfolding = 0;
        //start printing best switch of the group
        fprintf(textfile, "\n\nSwitch with highest mfe difference in group number %d\nLocation: %d\nSequence: ", (window-1), templocation);
        for (k = 0; k < strlen(bestswitch); k++){//print final switch sequence in txt file
            fprintf(textfile, "%c", bestswitch[k]);
        }
        makedna(bestswitch);
        fprintf(textfile, "\nDNA Sequence: ");
        for (k = 0; k < strlen(bestswitch); k++){//print final switch sequence in txt file
            fprintf(textfile, "%c", bestswitch[k]);
        }
        fprintf(textfile, "\nMFE: %fkcal/mol\nMFE difference: %fkcal/mol\n\n\n\n\n\n\n", mfetemp, mfedifference);
        //end printing best switch of the
        mfedifference = DBL_MIN;
    //
    fprintf(textfile, "\n\n\n\n\n\n\n\n\n\n\n\nHighest to lowest MFE difference\n\n\n\n");
    
    j = 0;
    
    for (x=0; x<switchnumber; x ++) {
        for (k=0; k<switchnumber; k++) {
            if (mfedifference<mfearray[k]) {
                mfedifference=mfearray[k];
                j = k;
            }
        }
        fprintf(textfile, "\nRank: %d\nSwitch Number: %d\nMFE difference: %fkcal/mol\n\n", x, j+1, mfedifference);
        mfedifference = DBL_MIN;
        mfearray[j] = DBL_MIN;
    }
    
    mfedifference = DBL_MIN;
    
    fprintf(textfile, "\n\n\n\n\n\n\n\n\n\n\n\nHighest to lowest MFE\n\n\n\n");
    
    j = 0;
    
    for (x=0; x<switchnumber; x ++) {
        for (k=0; k<switchnumber; k++) {
            printf("\nmfearray2: %f\n", mfearray2[k]);
            if (mfedifference<mfearray2[k]) {
                mfedifference=mfearray2[k];
                j = k;
            }
        }
        fprintf(textfile, "\nRank: %d\nSwitch Number: %d\nMFE difference: %fkcal/mol\n\n", x, j+1, mfedifference);
        mfedifference = DBL_MIN;
        mfearray2[j] = DBL_MIN;
    }
    

                                return 0;
}



int main(void) {

    int window = 0;
    double mfearray[len] = {'0'};
    char filename[len] = {'0'};
    char filenamecsv[len] = {'0'};
    char input[len] = {'0'};
    char promoterseq[100] = {'0'};
    char rbs[500] = {'0'};
    char res5[100] = {'0'};
    char res3[100] = {'0'};
    int repeat = 0;
    char copy[1200] = {'0'};
    int i, j, l;
    int n = 0;
    int between;
    int triggerlength = 0;
    int choice;//RBS and AUG choice
    int printstructure = 0;
    int plotrnachoice = 0;
    int promoter = 0;//T7 or not
    int maxdomain = 0;//maximum toehold domain base pairs
    int maxrbs = 0;//maximum rbs domain base pairs
    double inputtemperature = 0;//input temperature
    FILE *textfile;
    FILE *csvfile;
    
    
    vrna_md_defaults_temperature(VRNA_MODEL_DEFAULT_TEMPERATURE);//which is 37.0 degrees celsius

    
    //create data.txt file then writes output in
    printf("Hello\n");
    printf("Enter File Name: ");
    fgets(filename, 1000, stdin);
    strtok(filename, "\n");//removes '\n'
    textfile = fopen(filename, "w");
    strcpy(filenamecsv, filename);
    strcat(filenamecsv, " EXCEL.csv");
    csvfile = fopen(filenamecsv, "w");
    fprintf(textfile, "Hello!\n\n\n\n");
    
    do {
        printf("Enter DNA/RNA sequence(1000bp max at a time): ");
        fgets(input, len, stdin);
        strtok(input, "\n");//removes '\n' from input
        makerna(input);
    }while (check_valid(input) != 0);
    
    do {
        printf("Enter DNA/RNA sequence(1000bp max at a time)(Type enter if finished entering sequence): ");
        memset(copy,'\0',strlen(copy));
        fgets(copy, len, stdin);
        makerna(copy);
        strcat(input, copy);
        strtok(input, "\n");//removes '\n' from input
    }while (!(copy[0] == '\n'));
    
    memset(copy,'\0',strlen(copy));
    
    do {
        printf("Enter Promoter sequence(Type 'enter' for standard promoter(J23100){TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC})(Type '0' for T7 promoter{TAATACGACTCACTATA}): ");
        fgets(promoterseq, 99, stdin);
        if (promoterseq[0] == '\n') {
            strcpy(promoterseq, "TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC");
        }
        else if (promoterseq[0] == '0') {
            strcpy(promoterseq, "TAATACGACTCACTATA");
            promoter = 1;
            
        }
        strtok(promoterseq, "\n");//removes '\n' from input
        makerna(promoterseq);
    }while (check_valid(promoterseq) != 0);
    do {
        printf("Enter 5 end restriction site(Type 'enter' for standard 5 end restriction site(NotI+XbaI){GCGGCCGCUUCUAGA}): ");
        fgets(res5, 99, stdin);
        if (res5[0] == '\n') {
            strcpy(res5, "GCGGCCGCUUCUAGA");
        }
        makerna(res5);
        strtok(res5, "\n");//removes '\n' from rbs
        } while (check_valid(res5) != 0);
    do {
        printf("Enter 3 end restriction site(Type 'enter' for standard 3 end restriction site(XhoI){CUCGAG}): ");
        fgets(res3, 99, stdin);
        if (res3[0] == '\n') {
            strcpy(res3, "CUCGAG");
        }
        makerna(res3);
        strtok(res3, "\n");//removes '\n' from rbs
    } while (check_valid(res3) != 0);
    do {
        printf("Enter RBS sequence(Type 'enter' for standard rbs{AAAGAGGAGAAA}): ");
        fgets(rbs, 499, stdin);
        if (rbs[0] == '\n') {
            strcpy(rbs, "AAAGAGGAGAAA");
        }
        makerna(rbs);
        strtok(rbs, "\n");//removes '\n' from rbs
    } while (check_valid(rbs) != 0);
    do {
        //printf("Enter distance between groups: ");
        //scanf("%d", &between);
        between = 0;
        printf("Enter trigger length: ");
        scanf("%d", &triggerlength);
        printf("Enter maximum toehold domain base pairs: ");
        scanf("%d", &maxdomain);
        //printf("Enter maximum RBS domain base pairs: ");
        //scanf("%d", &maxrbs);
        maxrbs = INT_MAX;
        printf("Enter temperature in degrees celsius: ");
        scanf("%f", &inputtemperature);
        printf("Print MFE Structures? ('0' if not, '1' if yes): ");
        scanf("%d", &printstructure);
        printf("Plot RNA Picture? ('0' if not, '1' if yes): ");
        scanf("%d", &plotrnachoice);
        //printf("Eliminate for RBS and/or AUG pairings? ('0' for nocheck, '1' for RBS, '2' for AUG,'3' for BOTH): ");
        //scanf("%d", &choice);
        choice = 0;
        //printf("\n get\n");
    } while (choice_valid(choice) != 0);
    

    l = strlen(input);
    printf("\nstrlen: %d\nInput: %s\n\n", l, input);

    
    //printf("123\n");

    look(filename, mfearray, promoterseq, triggerlength, inputtemperature, maxdomain, maxrbs, promoter, choice, printstructure, plotrnachoice, between, res5, res3, input, rbs, l, window, n, textfile, csvfile);

    //printf("%d\n",round);
    


    fclose(textfile);
    fclose(csvfile);
    
    vrna_md_defaults_temperature(VRNA_MODEL_DEFAULT_TEMPERATURE);//which is 37.0 degrees celsius
    
    printf("-----------------------------------------\n");
    return 0;
}
