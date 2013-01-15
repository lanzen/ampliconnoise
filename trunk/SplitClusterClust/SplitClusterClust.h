typedef struct s_Params
{
  char* szDatFile;

  char* szUCFile;

  char* szMapFile;

  int nSplit;

  int nMinSize;
} t_Params;


typedef struct s_Map
{
  int nU;

  int *anNM;

  int **aanMap;

} t_Map;

typedef struct s_Data
{
  char **aszLines;
  
  int nN;

  int nM;
  
} t_Data;

/* Input Definitions */
#define OPTION  0      /* optional */
#define ALWAYS  1      /* required */

#define UC_FILE                 "-uin"
#define MAP_FILE                "-min"
#define DAT_FILE                "-din"
#define MIN_SIZE                "-m"

#define LIST_SUFFIX             ".list"
#define TREE_SUFFIX             ".tree"
#define DAT_SUFFIX              ".dat"

#define INTERNAL        -1
#define MAX_WORD_LENGTH 1024
#define MAX_LINE_LENGTH 1000000
#define MAX_SPLIT       10000
#define DELIM           " "
#define DELIMT          "\t"
#define BIG_LINE_LENGTH 400000000
#define DELIM2          ",:"

void getCommandLineParams(t_Params *ptParams,int argc,char *argv[]);

void readData(t_Data *ptData, char* szDatFile);

void writeData(FILE* ofp, t_Data *ptData, int nEntries, int *anEntries, t_Map *ptMap);

void readMapFile(t_Map *ptMap, char *szMapFile);

void readUCFile(int ***paanSplit, int *pnSplit, int **panN, char* szUCFile);
