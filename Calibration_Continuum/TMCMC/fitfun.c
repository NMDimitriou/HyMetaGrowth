#include "fitfun.h"


static pthread_mutex_t fork_mutex = PTHREAD_MUTEX_INITIALIZER;
static int flag[4096];  // MAX_WORKERS

#define REMOVEDIRS  0
#define PREFIX      "."  // "/scratch"
#define FAIL        -1e12
#define BUFLEN      1024




void fitfun_initialize() {
}



double fitfun(double *x, int n, void *output, int *winfo) {
    char workdir[BUFLEN], bindir[BUFLEN];
    double t, loglike;
    int i;

    // make tmp directory
    char cwd[BUFLEN]; getcwd(cwd, BUFLEN);
    sprintf(bindir, "%s/model", cwd);
    sprintf(workdir, PREFIX"/tmpdir.%d.%d.%d.%d",
            winfo[0], winfo[1], winfo[2], winfo[3]);
    mkdir(workdir, S_IRWXU);

    // info
    int pid = getpid();
    printf("Spawner(%d): running in %s with params", pid, workdir);
    for (i = 0; i < n; i++) printf(" %.6lf", x[i]);
    printf("\n");
    fflush(0);

    // fork simulation
    t = my_gettime();
    while (pthread_mutex_trylock(&fork_mutex) == EBUSY) usleep(500000);

    int rf = fork();
    if (rf < 0) {
        printf("Fork failed\n"); fflush(0);
    } else if (rf == 0) {



		chdir(workdir);

        // copy necessary stuff
        if (copy_from_dir(bindir) != 0) {
            printf("Error in copy from dir %s\n", bindir);
            abort();
        }


        // write parametes to the simulation's input file
        FILE *finp = fopen("params.txt", "w");
        for (i = 0; i < n; i++) fprintf(finp, "%.16lf\n", x[i]);
	fprintf(finp, "%d\n", winfo[1]);
        fprintf(finp, "%d"  , pid);
        fclose(finp);

        // run simulation
        char line[BUFLEN], *largv[64];
        sprintf(line, "sh ./doall.sh");
        parse(line, largv);
        int fd = open("output.txt", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
        dup2(fd, 1); dup2(fd, 2);  // make stdout and stderr go to file
        close(fd);

		execvp(*largv, largv);

		printf("This point must not be reached!\n");
		exit(1);
    }

	pthread_mutex_unlock(&fork_mutex);


    // wait for fork
    int status; waitpid(rf, &status, 0); sync();

    // read results
    char llfile[BUFLEN]; sprintf(llfile, "%s/loglike.txt", workdir);
    FILE * pFile = fopen(llfile, "r");
    if (pFile == NULL) {
        loglike = FAIL;
    } else {
        while (!feof(pFile)) fscanf(pFile, "%lf", &loglike);
        fclose (pFile);
    }
    if (isnan(loglike) || isinf(loglike)) loglike = FAIL;

    char mmfile[BUFLEN]; sprintf(mmfile, "%s/status.txt", workdir);
    char mstat[BUFLEN];
    const char* message="No status.";
    FILE * mFile = fopen(mmfile, "r");
    if (mFile == NULL) {
        sprintf(mstat,"%s", message);
    } else {
        while (!feof(mFile)) fgets(mstat, BUFLEN, mFile);
        fclose (mFile);
    }

    // cleanup
    if (REMOVEDIRS && flag[torc_worker_id()] == 0) rmrf(workdir);

    // info
    printf("task(%d.%d.%d.%d): ", winfo[0], winfo[1], winfo[2], winfo[3]);
    for (i = 0; i < n; i++) printf("%.6lf ", x[i]);
    printf(" = %.6lf in %lf secs", loglike, my_gettime()-t);
    printf(",%s\n",mstat);
    fflush(0);



    return loglike;
}


void fitfun_finalize() {
}
