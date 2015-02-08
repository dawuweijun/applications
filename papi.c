/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */

#include <papi.h>

	int EventSet = PAPI_NULL;

	PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_thread_init( (unsigned long (*)(void)) (omp_get_thread_num) );

    /* for master thread */
    PAPI_create_eventset(&EventSet);
    PAPI_add_event(EventSet, eventList[atoi(argv[2])]);
    PAPI_start(EventSet);
    PAPI_reset(EventSet);

    /* for other threads */
    #pragma omp parallel
    {

    	char string[1001];
		PAPI_event_code_to_name(eventList[atoi(argv[2])], string);
		int thread = omp_get_thread_num();

    	int MB_EventSet = PAPI_NULL;
	    PAPI_create_eventset(&MB_EventSet);
	    PAPI_add_event(MB_EventSet, eventList[atoi(argv[2])]);
		PAPI_start(MB_EventSet);

		#pragma omp for
		for

		PAPI_read(MB_EventSet, &valor);
        if(thread != 0)
			printf("%s => %d - %llu\n", string, thread, valor);
    }
    
    /* for master thread */
    char string[1001];
	PAPI_event_code_to_name(eventList[atoi(argv[2])], string);

	PAPI_read(EventSet, &valor);
	printf("%s => %d - %llu\n", string, 0, valor);
	PAPI_stop(EventSet, NULL);


	float rtime, ptime, mflops;
	long long flpops;
	PAPI_flops(&rtime, &ptime, &flpops, &mflops);
	printf( "\t%d - real time:       %f\n",thread, rtime);
	printf( "\t%d - process time:    %f\n",thread, ptime);
	printf( "\t%d - FP Operations:   %lld\n",thread, flpops);
	printf( "\t%d - MFLOPS           %f\n", thread, mflops);


	float rtime, ptime, ipc;
	long long ins;
	PAPI_ipc(&rtime, &ptime, &ins, &ipc);
	printf( "%d - real time:       %f\n",thread, rtime);
	printf( "%d - process time:    %f\n",thread, ptime);
	printf( "%d - Instructions:    %lld\n",thread, ins);
	printf( "%d - IPC              %f\n\n", thread, ipc);



	int PAPI_event_code_to_name(int  EventCode, char * EventName );
	int PAPI_event_name_to_code(char * EventName, int * EventCode );

	const int eventList[] = {PAPI_BR_INS};
	const int eventList[] = {PAPI_FAD_INS, PAPI_FDV_INS, PAPI_FMA_INS, PAPI_FML_INS, PAPI_FNV_INS, PAPI_FP_INS, PAPI_FP_OPS, PAPI_FSQ_INS};
	const int eventList[] = {PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_VEC_INS};
	const int eventList[] = {PAPI_L1_DCA, PAPI_L1_DCH, PAPI_L1_DCM, PAPI_L1_ICA, PAPI_L1_ICH, PAPI_L1_ICM, PAPI_L1_TCA, PAPI_L1_TCH, PAPI_L1_TCM};
	const int eventList[] = {PAPI_L2_DCA, PAPI_L2_DCH, PAPI_L2_DCM, PAPI_L2_ICA, PAPI_L2_ICH, PAPI_L2_ICM, PAPI_L2_TCA, PAPI_L2_TCH, PAPI_L2_TCM};
	const int eventList[] = {PAPI_L3_DCA, PAPI_L3_DCH, PAPI_L3_DCM, PAPI_L3_ICA, PAPI_L3_ICH, PAPI_L3_ICM, PAPI_L3_TCA, PAPI_L3_TCH, PAPI_L3_TCM};
	const int eventList[] = {PAPI_LD_INS, PAPI_LST_INS, PAPI_SR_INS};
	const int eventList[] = {PAPI_TLB_DM, PAPI_TLB_IM, PAPI_TLB_TL};