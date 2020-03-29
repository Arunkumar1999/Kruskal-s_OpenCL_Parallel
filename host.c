
#define PROGRAM_FILE "kernal.cl"
#pragma warning(disable: 4996)
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include<time.h>
#ifdef MAC
#include <OpenCL/cl.h>
#else  
#include <CL/cl.h>
#endif
typedef struct Edge
{
    int src, dest, weight;
}Edge;
void merge(struct Edge* arr, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    struct Edge* L = (struct Edge*)malloc(sizeof(struct Edge) * n1);
    struct Edge* R = (struct Edge*)malloc(sizeof(struct Edge) * n2);
         

    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 
    while (i < n1 && j < n2)
    {
        if (L[i].weight <= R[j].weight)
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort(struct Graph* arr, int l, int r)
{
    if (l < r)
    {
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        int m = l + (r - l) / 2;

        // Sort first and second halves 
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        merge(arr, l, m, r);
    }
}
// a structure to represent a connected, undirected 
// and weighted graph 
typedef struct Graph
{
    // V-> Number of vertices, E-> Number of edges 
    int V, E;

    // graph is represented as an array of edges. 
    // Since the graph is undirected, the edge 
    // from src to dest is also edge from dest 
    // to src. Both are counted as 1 edge here. 
    struct Edge* edge;
}Graph;

// Creates a graph with V vertices and E edges 
struct Graph* createGraph(int V, int E)
{
    struct Graph* graph = malloc(sizeof(Graph));
    graph->V = V;
    graph->E = E;

    graph->edge = malloc(sizeof(Edge) * E);

    return graph;
}

// A structure to represent a subset for union-find 
struct subset
{
    int parent;
    int rank;
};

// A utility function to find set of an element i 
// (uses path compression technique) 
int find(struct subset *subsets, int i)
{
    // find root and make root as parent of i 
    // (path compression) 
    //printf("%d in \n ", i);
    if (subsets[i].parent != i)
        subsets[i].parent = find(subsets, subsets[i].parent);

    return subsets[i].parent;
}

// A function that does union of two sets of x and y 
// (uses union by rank) 
void Union(struct subset subsets[], int x, int y)
{
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);

    // Attach smaller rank tree under root of high 
    // rank tree (Union by Rank) 
    if (subsets[xroot].rank < subsets[yroot].rank)
        subsets[xroot].parent = yroot;
    else if (subsets[xroot].rank > subsets[yroot].rank)
        subsets[yroot].parent = xroot;

    // If ranks are same, then make one as root and 
    // increment its rank by one 
    else
    {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}

// Compare two edges according to their weights. 
// Used in qsort() for sorting an array of edges 
int myComp(const void* a, const void* b)
{
    struct Edge* a1 = (struct Edge*)a;
    struct Edge* b1 = (struct Edge*)b;
    return a1->weight > b1->weight;
}

// The main function to construct MST using Kruskal's algorithm 
void KruskalMST(struct Graph* graph)
{
    int V = graph->V;
    struct Edge *result=malloc(sizeof(struct Edge)*(graph->E+1)); // Tnis will store the resultant MST 
    int e = 0; // An index variable, used for result[] 
    int i = 0; // An index variable, used for sorted edges 

    // Step 1: Sort all the edges in non-decreasing 
    // order of their weight. If we are not allowed to 
    // change the given graph, we can create a copy of 
    // array of edges 
    //qsort(graph->edge, graph->E, sizeof(graph->edge[0]), myComp);
    //mergeSort(graph->edge,0, graph->E-1);
    // Allocate memory for creating V ssubsets 
    struct subset* subsets =
        (struct subset*) malloc(V * sizeof(struct subset));
    //printf("%d\n", graph->E);
    //for (int i = 0; i < 65; i++)
    //{
      //  printf("%d %d %d %d iop\n", i, graph->edge[i].weight, graph->edge[i].src, graph->edge[i].dest);
    //}
    // Create V subsets with single elements 
    for (int v = 0; v < V; ++v)
    {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }

    // Number of edges to be taken is equal to V-1 
    while (e < V - 1 && i < graph->E)
    {
        //printf("%d %d %d %d weight\n",i, graph->edge[i].weight, graph->edge[i].src, graph->edge[i].dest);
        // Step 2: Pick the smallest edge. And increment 
        // the index for next iteration 
        struct Edge next_edge = graph->edge[i++];

        int x = find(subsets, next_edge.src);
        int y = find(subsets, next_edge.dest);

        // If including this edge does't cause cycle, 
        // include it in result and increment the index 
        // of result for next edge 
        if (x != y)
        {
            result[e++] = next_edge;
            Union(subsets, x, y);
        }
        // Else discard the next_edge 
    }

    // print the contents of result[] to display the 
    // built MST 
    FILE* ptr;
    ptr = fopen("sequential.txt", "w");
    
    printf("Following are the edges in the constructed MST\n");
    for (i = 0; i < e; ++i)
    {
        printf("%d -- %d == %d\n", result[i].src, result[i].dest,
            result[i].weight);
        fprintf(ptr, "%d %d %d\n", result[i].src, result[i].dest,
            result[i].weight);

    }
    fclose(ptr);
    return;
}

typedef struct node
{
    cl_uint src;
    cl_uint dest;
    cl_uint weight;
   
}node;
size_t getPaddedSize(size_t n)
{
    unsigned int log2val = (unsigned int)ceil(log((float)n) / log(2.f));
    return (size_t)pow(2, log2val);
}

size_t GetGlobalWorkSize(size_t DataElemCount, size_t LocalWorkSize)
{
    size_t r = DataElemCount % LocalWorkSize;
    if (r == 0)
        return DataElemCount;
    else
        return DataElemCount + LocalWorkSize - r;
}
void parallel_merge(struct node* arr, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    /* create temp arrays */
    struct node* L = (struct Edge*)malloc(sizeof(struct Edge) * n1);
    struct node* R = (struct Edge*)malloc(sizeof(struct Edge) * n2);


    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 
    while (i < n1 && j < n2)
    {
        if (L[i].weight <= R[j].weight)
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    /* Copy the remaining elements of L[], if there
       are any */
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }

    /* Copy the remaining elements of R[], if there
       are any */
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void parallel_mergeSort(struct node* arr, int l, int r)
{
    if (l < r)
    {
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        int m = l + (r - l) / 2;

        // Sort first and second halves 
        parallel_mergeSort(arr, l, m);
        parallel_mergeSort(arr, m + 1, r);

        parallel_merge(arr, l, m, r);
    }
}
int main() {

    /* Host/device data structures */
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;
    cl_int i, err;

    /* Program/kernel data structures */
    cl_program program;
    FILE* program_handle;
    char* program_buffer, * program_log;
    size_t program_size, log_size;
    cl_kernel kernel;

    /* Data and buffers */
    float mat[16], vec[4], result[4];
    float correct[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    cl_mem mat_buff, vec_buff, res_buff;
    size_t work_units_per_kernel;

    cl_mem            m_dPingArray;
    cl_mem            m_dPongArray;
    cl_mem            m_po;
    cl_mem            m_indicate;
    cl_mem            m_test;
    cl_mem            m_flag;
    cl_mem            m_size;
    cl_mem            sub_buffer;
    cl_kernel         m_BitonicStartKernel;
    cl_kernel         m_BitonicGlobalKernel;
    cl_kernel         m_BitonicLocalKernel;
    cl_kernel         kruskal_algo;
    cl_kernel         test_kruskal;
    size_t            m_N;
    size_t            m_N_padded;


    // input data
    //unsigned int* m_hInput;
    // results
    //unsigned int* m_resultCPU;
    //unsigned int* m_resultGPU[3];

    size_t LocalWorkSize[3] = { 256, 1, 1 };
    unsigned int arraySize = 4 * 4;

    m_N = arraySize;
    //m_N_padded = abs(getPaddedSize(m_N));
    m_N_padded = 1048576;
    printf("%d as\n", abs(m_N_padded));
    node* m_hInput = malloc(sizeof(node) * m_N_padded);
    node* m_resultCPU = malloc(sizeof(node) * m_N_padded);
    node* m_resultGPU[3];
    // input data
    //unsigned int m_hInput[m_N_padded];
    // results
    //unsigned int m_resultCPU[m_N_padded];
    //m_hInput = (unsigned int*)malloc(m_N_padded * sizeof(unsigned int));
    //m_resultCPU = (unsigned int*)malloc(m_N_padded * sizeof(unsigned int));
    //srand(time(0));
    //printf("%d asf %d %d \n", m_hInput == NULL,m_N,m_N_padded);
    //for (unsigned int i = 0; i < m_N_padded; i++)
    //{ //m_hInput[i] = m_N - i; 

    //    m_hInput[i].weight = rand();
    //}
    //printf("\n");
    int V = 1449; // Number of vertices in graph 
    int E = m_N_padded; // Number of edges in graph 
    struct Graph* sequential_input = createGraph(V,m_N_padded );

    int temp=0;
    for (cl_uint k = 0;k < 1450; k++)
    {
        for (cl_uint j = 0;j < k; j++)
        {   
            sequential_input->edge[temp].src = k;
            sequential_input->edge[temp].dest = j;

            m_hInput[temp].src = k;
            m_hInput[temp].dest = j;
            //printf("%d %d %d %d abc\n ", sequential_input->edge[temp].src, k,j, sequential_input->edge[temp].dest);
            m_hInput[temp].weight = rand() % 50+ rand() % 20;

           if (m_hInput[temp].weight > 55)m_hInput[temp].weight = 999;

           sequential_input->edge[temp].weight = m_hInput[temp].weight;
           temp++;
           //printf("%d \n", temp);
           if (temp == m_N_padded || temp== m_N_padded)
               break;
        }
        if (temp == m_N_padded)
            break;
        sequential_input->edge[temp].src = k;
        sequential_input->edge[temp].dest = k;
        sequential_input->edge[temp].weight = 999;

        m_hInput[temp].src = k;
        m_hInput[temp].dest = k;

        m_hInput[temp].weight = 999;
        temp++;
        if (temp == m_N_padded)
            break;
    }
    printf("%d %d hy\n", temp, m_hInput[0].src);
    m_N_padded = temp;
    //int temp = 0;
   // for (int k = 0;k < temp; k++)
    //{
    //    printf("%d %d %d input\n", sequential_input->edge[k].src, sequential_input->edge[k].dest, sequential_input->edge[k].weight);
    //}
    mergeSort(sequential_input->edge, 0, sequential_input->E - 1);
    parallel_mergeSort(m_hInput, 0, temp - 1);
    //m_hInput[0].src = 0;
    //m_hInput[0].dest = 1;
    //m_hInput[0].weight = 10;
    //m_hInput[0].indicate = 0;

    // add edge 0-2  
    //m_hInput[1].src = 0;
    //m_hInput[1].dest = 2;
    //m_hInput[1].weight = 6;
    //m_hInput[1].indicate = 0;

    // add edge 0-3  
    //m_hInput[2].src = 0;
    //m_hInput[2].dest = 3;
    //m_hInput[2].weight = 5;
    //m_hInput[2].indicate = 0;
    // add edge 1-3  
    //m_hInput[3].src = 1;
    //m_hInput[3].dest = 3;
    //m_hInput[3].weight = 15;
   // m_hInput[3].indicate = 0;
    // add edge 2-3  
    //m_hInput[4].src = 2;
    //m_hInput[4].dest = 3;
    //m_hInput[4].weight = 4;
  //  m_hInput[4].indicate = 0;
    double time_spent = 0.0;
    clock_t begin = clock();
    KruskalMST(sequential_input);
    clock_t end = clock();
    time_spent = time_spent + (double)(end - begin) / CLOCKS_PER_SEC;
    printf("sequrntial total time %f\n", time_spent);
    time_spent = 0.0;
     //begin = clock();

    printf("hkdvkh\n");
    /* Identify a platform */
    err = clGetPlatformIDs(1, &platform, NULL);
    if (err < 0) {
        perror("Couldn't find any platforms");
        exit(1);
    }

    /* Access a device */
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
    if (err < 0) {
        perror("Couldn't find any devices");
        exit(1);
    }

    /* Create the context */
    context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    if (err < 0) {
        perror("Couldn't create a context");
        exit(1);
    }

    /* Read program file and place content into buffer */
    program_handle = fopen(PROGRAM_FILE, "rb");
    if (program_handle == NULL) {
        perror("Couldn't find the program file");
        exit(1);
    }
    fseek(program_handle, 0, SEEK_END);
    program_size = ftell(program_handle);
    rewind(program_handle);
    program_buffer = (char*)malloc(program_size + 1);
    program_buffer[program_size] = '\0';
    fread(program_buffer, sizeof(char), program_size, program_handle);
    fclose(program_handle);

    /* Create program from file */
    program = clCreateProgramWithSource(context, 1,
        (const char**)&program_buffer, &program_size, &err);
    if (err < 0) {
        perror("Couldn't create the program");
        exit(1);
    }
    free(program_buffer);

    /* Build program */
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err < 0) {

        /* Find size of log and print to std output */
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
            0, NULL, &log_size);
        program_log = (char*)malloc(log_size + 1);
        program_log[log_size] = '\0';
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
            log_size + 1, program_log, NULL);
        printf("%s\n", program_log);
        free(program_log);
        exit(1);
    }

    cl_int clError, clError2;
    m_dPingArray = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(node) * m_N_padded, NULL, &clError2);
    clError = clError2;
    m_dPongArray = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(node) * m_N_padded, NULL, &clError2);
    clError |= clError2;
    m_po = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(node) * m_N_padded, NULL, &clError2);
    clError |= clError2;
    m_indicate = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_uint) * m_N_padded, NULL, &clError2);
    clError |= clError2;
    m_test = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_uint) * m_N_padded, NULL, &clError2);
    clError |= clError2;
    m_flag = clCreateBuffer(context, CL_MEM_READ_WRITE, 1 * sizeof(cl_uint), NULL, &clError2);
    clError |= clError2;
    m_size = clCreateBuffer(context, CL_MEM_READ_WRITE, 1 * sizeof(cl_uint), NULL, &clError2);
    clError |= clError2;

    if (clError < 0) {
        printf("pm %d", clError);
        exit(1);
    }

    m_BitonicStartKernel = clCreateKernel(program, "Sort_BitonicMergesortStart", &clError);
    if (clError < 0)
    {
        printf("opiuh %d", clError);
        exit(1);
    }
    m_BitonicGlobalKernel = clCreateKernel(program, "Sort_BitonicMergesortGlobal", &clError);
    if (clError < 0)
    {
        printf("opiu %d", clError);
        exit(1);
    }
    m_BitonicLocalKernel = clCreateKernel(program, "Sort_BitonicMergesortLocal", &clError);
    if (clError < 0) {
        printf("qweer %d", clError);
        exit(1);
    }
    kruskal_algo = clCreateKernel(program, "kruskal_algo", &clError);
    if (clError < 0) {
        printf("qweer %d", clError);
        exit(1);
    }
    test_kruskal = clCreateKernel(program, "test_kruskal", &clError);
    if (clError < 0) {
        printf("qweer %d", clError);
        exit(1);
    }

    size_t globalWorkSize[1];
    size_t localWorkSize[1];

    localWorkSize[0] = LocalWorkSize[0];
    globalWorkSize[0] = GetGlobalWorkSize(m_N_padded / 2, localWorkSize[0]);

    unsigned int limit = (unsigned int)2 * LocalWorkSize[0]; //limit is double the localWorkSize

    // start with Sort_BitonicMergesortLocalBegin to sort local until we reach the limit
    clError = clSetKernelArg(m_BitonicStartKernel, 0, sizeof(cl_mem), (void*)&m_dPingArray);
    clError |= clSetKernelArg(m_BitonicStartKernel, 1, sizeof(cl_mem), (void*)&m_dPongArray);

    if (clError < 0) {
        printf("op %d", clError);
        exit(1);
    }
    queue = clCreateCommandQueue(context, device, 0, &err);
    if (err < 0) {
        perror("Couldn't create the command queue");
        exit(1);
    }
    clEnqueueWriteBuffer(queue, m_dPingArray, CL_FALSE, 0, m_N_padded * sizeof(node), m_hInput, 0, NULL, NULL);

    clError = clEnqueueNDRangeKernel(queue, m_BitonicStartKernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (clError < 0)
    {
        perror("Couldn't enqueue the kernel execution command");
        exit(1);
    }

    for (unsigned int blocksize = limit; blocksize <= m_N_padded; blocksize <<= 1) {
        for (unsigned int stride = blocksize / 2; stride > 0; stride >>= 1) {
            if (stride >= limit) {
                //Sort_BitonicMergesortGlobal
                clError = clSetKernelArg(m_BitonicGlobalKernel, 0, sizeof(cl_mem), (void*)&m_dPongArray);
                clError |= clSetKernelArg(m_BitonicGlobalKernel, 1, sizeof(cl_uint), (void*)&m_N_padded);
                clError |= clSetKernelArg(m_BitonicGlobalKernel, 2, sizeof(cl_uint), (void*)&blocksize);
                clError |= clSetKernelArg(m_BitonicGlobalKernel, 3, sizeof(cl_uint), (void*)&stride);
                if (clError < 0) {
                    printf("was %d", clError);
                    exit(1);
                }

                clError = clEnqueueNDRangeKernel(queue, m_BitonicGlobalKernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
                if (clError < 0)
                    exit(1);
            }
            else {
                //Sort_BitonicMergesortLocal
                clError = clSetKernelArg(m_BitonicLocalKernel, 0, sizeof(cl_mem), (void*)&m_dPongArray);
                clError |= clSetKernelArg(m_BitonicLocalKernel, 1, sizeof(cl_uint), (void*)&m_N_padded);
                clError |= clSetKernelArg(m_BitonicLocalKernel, 2, sizeof(cl_uint), (void*)&blocksize);
                clError |= clSetKernelArg(m_BitonicLocalKernel, 3, sizeof(cl_uint), (void*)&stride);

                if (clError < 0)
                    exit(1);
                clError = clEnqueueNDRangeKernel(queue, m_BitonicLocalKernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
                if (clError < 0)
                {
                    printf("qw %d", clError);
                    exit(1);
                }
            }
        }
    }
    //printf("school");
    //swap(m_dPingArray, m_dPongArray);
    m_resultGPU[0] = (node*)malloc(m_N_padded * sizeof(node));
    clError = clEnqueueReadBuffer(queue, m_dPongArray, CL_TRUE, 0, m_N_padded * sizeof(node), m_resultGPU[0], 0, NULL, NULL);
    if (clError < 0)
    {
        printf("as %d", clError);
        exit(1);
    }
    m_resultGPU[1] = (node*)malloc(m_N_padded * sizeof(node));
    for (unsigned int i = 0; i <m_N_padded; i++) {
        m_resultGPU[1][m_N_padded-i-1].src=m_hInput[i].src;
        m_resultGPU[1][m_N_padded - i - 1].dest = m_hInput[i].dest;
        m_resultGPU[1][m_N_padded - i - 1].weight = m_hInput[i].weight;
        //printf("%d %d %d  %d we\n", m_resultGPU[0][i].src, m_resultGPU[0][i].dest, m_resultGPU[0][i].weight, i);
        //cout << m_resultGPU[Task][i] << ",\t";
    }
    cl_uint* p = malloc(sizeof(cl_uint) * m_N_padded);
    cl_uint* indicate = malloc(sizeof(cl_uint) * m_N_padded);
    cl_uint* result_indicate = malloc(sizeof(cl_uint) * m_N_padded);
    cl_uint* test = malloc(sizeof(cl_uint) * m_N_padded);
    cl_uint* flag = malloc(sizeof(cl_uint) * m_N_padded);
    node* inter_m = malloc(sizeof(node) * m_N_padded);
    inter_m[m_N_padded - 1].src = m_resultGPU[1][m_N_padded-1].src;
    inter_m[m_N_padded - 1].dest = m_resultGPU[1][m_N_padded - 1].dest;
    flag[0] = 0;
    flag[1] = 0;
    flag[2] = 1;
    flag[3] = 1;
    flag[4] = 1;
    for (int i = 0; i < m_N_padded; i++)
    {
        p[i] = i;
        indicate[i] = 0;
    }

    //clError = clSetKernelArg(kruskal_algo, 0, sizeof(cl_mem), (void*)&m_dPingArray);
    //clError |= clSetKernelArg(kruskal_algo, 1, sizeof(cl_mem), (void*)&m_dPongArray);
    //clError |= clSetKernelArg(kruskal_algo, 2, sizeof(cl_mem), (void*)&m_po);
    //clError |= clSetKernelArg(kruskal_algo, 3, sizeof(cl_mem), (void*)&m_indicate);
    //clError |= clSetKernelArg(kruskal_algo, 4, sizeof(cl_mem), (void*)&m_test);
    //clError |= clSetKernelArg(kruskal_algo, 5, sizeof(cl_mem), (void*)&m_flag);
    //if (clError < 0) {
     //   printf("op %d", clError);
       // exit(1);
    //}
    //clEnqueueWriteBuffer(queue, m_dPingArray, CL_FALSE, 0, m_N_padded * sizeof(node), m_resultGPU[0], 0, NULL, NULL);
    //clEnqueueWriteBuffer(queue, m_po, CL_FALSE, 0, 5 * sizeof(cl_uint),p, 0, NULL, NULL);
    //clEnqueueWriteBuffer(queue, m_indicate, CL_FALSE, 0, m_N_padded * sizeof(cl_uint), indicate, 0, NULL, NULL);
    //clEnqueueWriteBuffer(queue, m_flag, CL_FALSE, 0, 5*sizeof(cl_uint), flag, 0, NULL, NULL);
    globalWorkSize[0] = m_N_padded;

    //clError = clEnqueueNDRangeKernel(queue,kruskal_algo, 1, NULL, globalWorkSize, globalWorkSize, 0, NULL, NULL);
    //if (clError < 0)
    //{
     //   perror("Couldn't enqueue the kernel execution command");
      //  exit(1);
    //}
    //m_resultGPU[0] = (node*)malloc(m_N_padded * sizeof(node));
    //clError = clEnqueueReadBuffer(queue, m_dPongArray, CL_TRUE, 0, m_N_padded* sizeof(node), m_resultGPU[0], 0, NULL, NULL);
    //clError = clEnqueueReadBuffer(queue, m_indicate, CL_TRUE, 0, m_N_padded * sizeof(cl_uint), result_indicate, 0, NULL, NULL);
    //clError = clEnqueueReadBuffer(queue, m_test, CL_TRUE, 0, m_N_padded * sizeof(cl_uint), test, 0, NULL, NULL);
    //clError = clEnqueueReadBuffer(queue, m_po, CL_TRUE, 0, m_N_padded * sizeof(cl_uint), p, 0, NULL, NULL);
    //clError = clEnqueueReadBuffer(queue, m_flag, CL_TRUE, 0, 5 * sizeof(cl_uint), flag, 0, NULL, NULL);

    //if (clError < 0)
    //{
      //  printf("as %d", clError);
       // exit(1);
    //}
    //for (unsigned int i = 0; i < 2; i++) {
      //  printf("%d %d %d %d %d %d fgfhgf\n", m_resultGPU[0][i].src, m_resultGPU[0][i].dest, result_indicate[i],test[i],p[i], i);

    //}
    //for (unsigned int i = 0; i < 5; i++) {
      //  printf("%d anbs\n",flag[i]);

    //}
    cl_mem m_temp_input;
    m_temp_input = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(node) * m_N_padded, m_resultGPU[1], &clError2);
    //m_dPingArray = clCreateBuffer(context, CL_MEM_READ_WRITE| CL_MEM_COPY_HOST_PTR, sizeof(node) * m_N_padded, NULL, &clError2);
    clError = clSetKernelArg(test_kruskal, 0, sizeof(cl_mem), (void*)&m_temp_input);
    
    clError |= clSetKernelArg(test_kruskal, 1, sizeof(cl_mem), (void*)&m_po);
    clError |= clSetKernelArg(test_kruskal, 2, sizeof(cl_mem), (void*)&m_indicate);
    clError |= clSetKernelArg(test_kruskal, 3, sizeof(cl_mem), (void*)&m_test);
    clError |= clSetKernelArg(test_kruskal, 4, sizeof(cl_mem), (void*)&m_size);
    clError |= clSetKernelArg(test_kruskal, 5, sizeof(cl_mem), (void*)&m_dPongArray);
    if (clError < 0) {
        printf("op %d", clError);
        exit(1);
    }
    cl_buffer_region region;


    begin = clock();
    int v, j, ml;
    globalWorkSize[0] = m_N_padded;
    //cl_buffer_region region;
    for (int o = m_N_padded-1;o >=0 ;o--)
    //for (int o = 0;o < m_N_padded ;o++)
    {
      // printf("%d %d %d %d ert \n", p[2],p[0], m_resultGPU[0][o].weight,indicate[5]);
        if (m_resultGPU[1][o].weight != 999 && indicate[o] != 2)
        {
            //printf("%d  ert \n", indicate[o]);
            //v = m_resultGPU[0][o].src;

            //while (p[v] != v)
            //{
              //  v = p[v];
            //}

            ml = inter_m[o].src;
            //printf("%d %d %d %dinter\n", ml, inter_m[o].src, inter_m[o].dest);
            //printf("%d op", o);
            //v = m_resultGPU[0][o].dest;
            //while (p[v] != v)
            //{
              //  v = p[v];
            //}
            j = inter_m[o].dest;
            //printf("%d %d %d %dinter\n", ml, inter_m[o].src,j, inter_m[o].dest);
            if (ml != j)
            {
                //printf("%d %d %d %d %d rtr\n", ml, j, m_resultGPU[0][o].src, m_resultGPU[0][o].dest,o);
                indicate[o] = 1;
                if (j > ml)
                    p[j] = ml;
                else
                    p[ml] = j;
            }
            else
            {
                indicate[o] = 2;
            }


            if (o>0)
            //if (o != m_N_padded - 1)
            {

                //region.origin = (o + 1) * sizeof(node);
                //region.size = (m_N_padded - o - 1) * sizeof(node);
               //sub_buffer = clCreateSubBuffer(m_dPingArray, CL_MEM_READ_WRITE, CL_BUFFER_CREATE_TYPE_REGION, &region, &clError);
                //clEnqueueWriteBuffer(queue, m_dPingArray, CL_FALSE, 0, o * sizeof(node), m_resultGPU[1], 0, NULL, NULL);
                //clEnqueueWriteBuffer(queue, m_dPingArray, CL_FALSE, 0, m_N_padded * sizeof(node), m_resultGPU[0], 0, NULL, NULL);
                clEnqueueWriteBuffer(queue, m_po, CL_FALSE, 0, m_N_padded * sizeof(cl_uint), p, 0, NULL, NULL);
                clEnqueueWriteBuffer(queue, m_indicate, CL_FALSE, 0, m_N_padded * sizeof(cl_uint), indicate, 0, NULL, NULL);
                clEnqueueWriteBuffer(queue, m_size, CL_FALSE, 0, sizeof(cl_uint), &o, 0, NULL, NULL);
                //globalWorkSize[0] = m_N_padded;
                globalWorkSize[0] = o;
                //localWorkSize[0] = LocalWorkSize[0];

                clError = clEnqueueNDRangeKernel(queue, test_kruskal, 1, NULL, globalWorkSize, NULL, 0, NULL, NULL);
                if (clError < 0)
                {
                    perror("Couldn't enqueue the kernel execution command");
                    exit(1);
                }
                clError = clEnqueueReadBuffer(queue, m_indicate, CL_TRUE, 0, m_N_padded * sizeof(cl_uint), indicate, 0, NULL, NULL);
                clError = clEnqueueReadBuffer(queue, m_dPongArray, CL_TRUE, 0, o * sizeof(node), inter_m, 0, NULL, NULL);
                //clError = clEnqueueReadBuffer(queue, m_dPongArray, CL_TRUE, 0, m_N_padded * sizeof(cl_uint), inter_m, 0, NULL, NULL);
                //clError = clEnqueueReadBuffer(queue, m_test, CL_TRUE, 0, 5 * sizeof(cl_uint), test, 0, NULL, NULL);
                //printf("%d %d %d %d %d %d wer\n", test[0], test[1],p[2],p[1],indicate[5],o);
                if (clError < 0)
                {
                    perror("Couldn't enqueue the kernel execution command");
                    exit(1);
                }
                //for (unsigned int z = 0; z < m_N_padded; z++) {
                //    printf("%d %d qert \n", indicate[z],o);
               // }
            }
        }
    }
     end = clock();
     FILE* ptr;
     ptr = fopen("parallel.txt", "w");
    for ( int u = m_N_padded-1; u >=0; u--) {
        //printf("%d \n ", u);
        if(indicate[u]==1){
        printf(" %d -->   %d  , weight=%d \n", m_resultGPU[1][u].src, m_resultGPU[1][u].dest, m_resultGPU[1][u].weight);
        fprintf(ptr, "%d %d %d\n", m_resultGPU[1][u].src, m_resultGPU[1][u].dest, m_resultGPU[1][u].weight);
        }
    
    
    }
    fclose(ptr);
    time_spent = time_spent + (double)(end - begin) / CLOCKS_PER_SEC;
    printf("parallel total time %f for %d size edge \n", time_spent,m_N_padded);
    clReleaseMemObject(m_dPingArray);
    clReleaseMemObject(m_dPongArray);
    clReleaseMemObject(m_po);
    clReleaseMemObject(m_indicate);
    clReleaseMemObject(m_test);
    clReleaseKernel(test_kruskal);
    clReleaseKernel(m_BitonicGlobalKernel);
    clReleaseKernel(m_BitonicLocalKernel);
    clReleaseKernel(m_BitonicStartKernel);
    clReleaseProgram(program);
    clReleaseContext(context);
    clReleaseCommandQueue(queue);
}



