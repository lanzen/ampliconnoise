#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#include "SplitClusterClust.h"

static char *usage[] = {"SplitClusterClust\n",
			"-din    string  dat filename\n",
			"-min    string  map filename\n",
			"-uin    string  uclust filename\n",
			"-m      integer min size\n"};

static int nLines = 5;

int main(int argc, char* argv[])
{
  FILE      *ifp    = NULL;
  t_Params  tParams;
  t_Data    tData;
  t_Map     tMap;
  int       **aanSplit = NULL;
  int        nSplit = 0;
  int       *anN = NULL;
  int        iL = 0, nLast = 0, nCount = 0, i = 0, j = 0;
  int        *anLast = NULL;
  char       szDir[MAX_WORD_LENGTH];
  char       szDatFile[MAX_WORD_LENGTH];
  FILE       *dfp = NULL;

  getCommandLineParams(&tParams, argc, argv);

  readData(&tData, tParams.szDatFile);

  readMapFile(&tMap, tParams.szMapFile);

  readUCFile(&aanSplit, &nSplit, &anN, tParams.szUCFile);

  for(i = 0; i < nSplit; i++){
    if(anN[i] < tParams.nMinSize){
      nLast += anN[i];
    }
  }

  i = 0;
  for(i = 0; i < nSplit; i++){ 
    if(anN[i] >= tParams.nMinSize){
      sprintf(szDir, "C%03d",nCount);
    
      mkdir(szDir, S_IRWXU);

      sprintf(szDatFile,"%s/%s%s",szDir,szDir,DAT_SUFFIX);
    
      printf("%d %d\n",i,anN[i]);

      dfp = fopen(szDatFile, "w");

      if(dfp){
	writeData(dfp, &tData, anN[i], aanSplit[i], &tMap);
    
	fclose(dfp);
      }

      nCount++;
    }
  }
  
  if(nLast > 0){
    anLast = (int *) malloc(sizeof(int)*nLast);
    
    printf("%d %d\n",i,nLast);
    iL = nCount;
    nCount=0;
    for(i = 0; i < nSplit; i++){
      if(anN[i] < tParams.nMinSize){
	for(j = 0; j < anN[i]; j++){
	  anLast[nCount + j] = aanSplit[i][j];
	}
	nCount += anN[i];
      }
    }

    if(nCount > 0){
      sprintf(szDir, "C%03d+",iL);
  
      mkdir(szDir, S_IRWXU);

      sprintf(szDatFile,"%s/%s%s",szDir,szDir,DAT_SUFFIX);

      dfp = fopen(szDatFile, "w");

      if(dfp){
	writeData(dfp, &tData, nLast, anLast, &tMap);
    
	fclose(dfp);
      }
    }

    free(anLast);
  }
  exit(EXIT_SUCCESS);
}

void writeUsage(FILE* ofp)
{
  int i = 0;
  char *line;

  for(i = 0; i < nLines; i++){
    line = usage[i];
    fputs(line,ofp);
  }
}

char *extractParameter(int argc, char **argv, char *param,int when)
{
  int i = 0;

  while((i < argc) && (strcmp(param,argv[i]))){
    i++;
  }

  if(i < argc - 1){
    return(argv[i + 1]);
  }

  if((i == argc - 1) && (when == OPTION)){
    return "";
  }

  if(when == ALWAYS){
    fprintf(stdout,"Can't find asked option %s\n",param);
  }

  return (char *) NULL;
}

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[])
{
  char *szTemp = NULL;
  char *pcError = NULL;

  /*get parameter file name*/
  ptParams->szDatFile  = extractParameter(argc,argv,DAT_FILE,ALWAYS);  
  if(ptParams->szDatFile == NULL)
    goto error;

  ptParams->szUCFile  = extractParameter(argc,argv,UC_FILE,ALWAYS);  
  if(ptParams->szUCFile == NULL)
    goto error;

  ptParams->szMapFile  = extractParameter(argc,argv,MAP_FILE, ALWAYS);  
  if(ptParams->szMapFile == NULL)
    goto error;

  szTemp = extractParameter(argc,argv,MIN_SIZE,ALWAYS);
  if(szTemp != NULL){
    ptParams->nMinSize = strtol(szTemp,&pcError,10);
    if(*pcError != '\0'){
      goto error;
    }
  }
  else{
    goto error;
  }

  return;

 error:
  writeUsage(stdout);
  exit(EXIT_FAILURE);

}




void readData(t_Data *ptData, char* szDatFile)
{
  char szLine[MAX_LINE_LENGTH];
  char *szTok = NULL, *pcError = NULL, *szBrk = NULL;;
  FILE *ifp = NULL;
  int i = 0, nM = 0, nN = 0;

  ifp = fopen(szDatFile, "r");
  
  if(ifp){
    fgets(szLine, MAX_LINE_LENGTH, ifp);

    //szBrk = strpbrk(szLine, "\n");

    szTok = strtok(szLine, DELIM);

    nN = strtol(szTok, &pcError, 10);
    if(*pcError != '\0')
      goto fileFormatError;

    szTok = strtok(NULL, DELIM);

    nM = strtol(szTok, &pcError, 10);
    if(*pcError != '\n')
      goto fileFormatError;

    ptData->nN = nN; ptData->nM = nM;

    ptData->aszLines = (char **) malloc(sizeof(char *)*3*nN);

    for(i = 0; i < 3*nN; i++){
      
      fgets(szLine, MAX_LINE_LENGTH, ifp);

      ptData->aszLines[i] = strdup(szLine);
    }

    fclose(ifp);
  }
  else{
    printf("Failed to open %s\n",szDatFile);
  }
  
  return;

 fileFormatError:
  fprintf(stderr, "Incorrectly formatted input file in readMeans");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void writeData(FILE* ofp, t_Data *ptData, int nEntries, int *anEntries, t_Map *ptMap)
{
  int i = 0, j = 0, index = 0;
  int nTotal = 0;

  for(i = 0; i < nEntries; i++){
    index = anEntries[i];
    nTotal += ptMap->anNM[index];
  }

  fprintf(ofp, "%d %d\n",nTotal,ptData->nM);

  for(i = 0; i < nEntries; i++){
    int index = anEntries[i];

    for(j = 0; j < ptMap->anNM[index]; j++){
      int nI = ptMap->aanMap[index][j];

      //fprintf(ofp, "%d ", nI);
      fputs (ptData->aszLines[nI], ofp);
    }
  }
}



void readMapFile(t_Map *ptMap, char *szMapFile)
{
  int i = 0, j = 0;
  char *szLine  = (char *) malloc(BIG_LINE_LENGTH*sizeof(char));
  FILE *ifp     = NULL;
  char *pcError = NULL, *szTemp  = NULL, *szTok = NULL;


  ifp = fopen(szMapFile, "r");
  
  if(ifp){
    
    ptMap->nU = 0;
    while(fgets(szLine, BIG_LINE_LENGTH, ifp) != NULL){
      ptMap->nU++;
    }
    fclose(ifp);

    ifp = fopen(szMapFile, "r");

    ptMap->anNM  = (int *) malloc(ptMap->nU*sizeof(int));
    ptMap->aanMap = (int **) malloc(ptMap->nU*sizeof(int *));
    
    if(!ptMap->anNM || !ptMap->aanMap)
      goto memoryError;

    for(i = 0; i < ptMap->nU; i++){
      fgets(szLine, BIG_LINE_LENGTH, ifp);

      szTok = strtok(szLine, DELIM2);

      for(j = 0; j < 3; j++){
	szTok = strtok(NULL, DELIM2);
      }

      ptMap->anNM[i] = strtol(szTok,&pcError, 10);
     
      ptMap->aanMap[i] = (int *) malloc(ptMap->anNM[i]*sizeof(int));
      if(!ptMap->aanMap[i])
	goto memoryError;

      for(j = 0; j < ptMap->anNM[i]; j++){
	szTok = strtok(NULL, DELIM2);
	
	ptMap->aanMap[i][j] = strtol(szTok,&pcError, 10);
      }
      
    }


    fclose(ifp);
  }
  else{
    goto fileError;
  }
  
  free(szLine);

  return;

 fileError:
  fprintf(stderr, "Failed to open map file\n", szMapFile);
  exit(EXIT_FAILURE);
 memoryError:
  fprintf(stderr, "Failed to allocate memory\n");
  exit(EXIT_FAILURE);
}

void readUCFile(int ***paanSplit, int *pnSplit, int **panN, char* szUCFile)
{
  char szLine[MAX_LINE_LENGTH];
  char *szTok = NULL, *pcError = NULL, *szBrk = NULL;;
  FILE *ifp = NULL;
  int i = 0, j = 0;
  int **aanSplit = NULL;
  int nSplit = 0;
  int *anN = (int *) malloc(sizeof(int)*MAX_SPLIT);
  int nID = -1;
  for(i = 0; i < MAX_SPLIT; i++){
    anN[i] = 0;
  }

  if(!anN)
    goto memoryError;

  ifp = fopen(szUCFile, "r");
  
  if(ifp){
    while(fgets(szLine, MAX_LINE_LENGTH, ifp)!=NULL){
      szTok = strtok(szLine, DELIMT);

      if(*szTok == 'S'){
	nSplit++;
      }

      if(*szTok == 'S' || *szTok == 'H'){
	int nI = -1;
	szTok = strtok(NULL, DELIMT);
	nI = strtol(szTok,&pcError,10);
	if(*pcError != '\0')
	  goto fileFormatError;
	
	anN[nI]++;
      }
    }

    fclose(ifp);
  }
  else{
    printf("Failed to open %s\n",szUCFile);
  }
  
  printf("%d\n",nSplit);

  aanSplit = (int **) malloc(nSplit*sizeof(int*));
  if(!aanSplit)
    goto memoryError;

  for(i = 0; i < nSplit; i++){
    aanSplit[i] = (int *) malloc(anN[i]*sizeof(int));
    anN[i] = 0;
    if(!aanSplit[i])
      goto memoryError;
  }

  ifp = fopen(szUCFile, "r");
  
  if(ifp){
    while(fgets(szLine, MAX_LINE_LENGTH, ifp)!=NULL){
      szTok = strtok(szLine, DELIMT);

      if(*szTok == 'S' || *szTok == 'H'){
	int nI = -1;
	szTok = strtok(NULL, DELIMT);

	nI = strtol(szTok,&pcError,10);
	if(*pcError != '\0')
	  goto fileFormatError;

	for(j = 0; j < 7; j++)
	  szTok = strtok(NULL, DELIMT);
	  
	nID = strtol(szTok,&pcError,10);
	if(*pcError != '\0')
	  goto fileFormatError;

	aanSplit[nI][anN[nI]] = nID;

	anN[nI]++;

      }
    }

    fclose(ifp);
  }
  else{
    printf("Failed to open %s\n",szUCFile);
  }

  (*paanSplit) = aanSplit;
  (*pnSplit) = nSplit;
  (*panN) = anN;
  return;

 fileFormatError:
  fprintf(stderr, "Incorrectly formatted input file in readUCFile");
  fflush(stderr);
  exit(EXIT_FAILURE);

 memoryError:
  fprintf(stderr, "Failed allocating memory in readUCFile");
  fflush(stderr);
  exit(EXIT_FAILURE);
}



