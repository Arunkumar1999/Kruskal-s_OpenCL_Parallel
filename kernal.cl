#define MAX_LOCAL_SIZE 256

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
typedef struct node
{
    unsigned int src;
    unsigned dest;
   unsigned weight;
   
}node;

inline void swap(node *a, node *b) {
	node tmp;
	tmp = *b;
	*b = *a;
	*a = tmp;
}

// dir == 1 means ascending
inline void sort(node *a, node *b, char dir) {
	if ((a->weight > b->weight) == dir) swap(a, b);
}

inline void swapLocal(__local node *a, __local node *b) {
	node tmp;
	tmp = *b;
	*b = *a;
	*a = tmp;
}

// dir == 1 means ascending
inline void sortLocal(__local node *a, __local node *b, char dir) {
	if ((a->weight > b->weight) == dir) swapLocal(a, b);
}

__kernel void Sort_BitonicMergesortStart(const __global node* inArray, __global node* outArray)
{
	__local node local_buffer[MAX_LOCAL_SIZE * 2];
	const uint gid = get_global_id(0);
	const uint lid = get_local_id(0);

	uint index = get_group_id(0) * (MAX_LOCAL_SIZE * 2) + lid;
	//load into local mem
	local_buffer[lid].src = inArray[index].src;
	local_buffer[lid].dest = inArray[index].dest;
	local_buffer[lid].weight = inArray[index].weight;

	local_buffer[lid + MAX_LOCAL_SIZE].src = inArray[index + MAX_LOCAL_SIZE].src;
	local_buffer[lid + MAX_LOCAL_SIZE].dest = inArray[index + MAX_LOCAL_SIZE].dest;
		local_buffer[lid + MAX_LOCAL_SIZE].weight = inArray[index + MAX_LOCAL_SIZE].weight;
	uint clampedGID = gid & (MAX_LOCAL_SIZE - 1);

	// bitonic merge
	for (uint blocksize = 2; blocksize < MAX_LOCAL_SIZE * 2; blocksize <<= 1) {
		char dir = (clampedGID & (blocksize / 2)) == 0; // sort every other block in the other direction (faster % calc)
#pragma unroll
		for (uint stride = blocksize >> 1; stride > 0; stride >>= 1){
			barrier(CLK_LOCAL_MEM_FENCE);
			uint idx = 2 * lid - (lid & (stride - 1)); //take every other input BUT starting neighbouring within one block
			sortLocal(&local_buffer[idx], &local_buffer[idx + stride], dir);
		}
	}

	// bitonic merge for biggest group is special (unrolling this so we dont need ifs in the part above)
	char dir = (clampedGID & 0); //even or odd? sort accordingly
#pragma unroll
	for (uint stride = MAX_LOCAL_SIZE; stride > 0; stride >>= 1){
		barrier(CLK_LOCAL_MEM_FENCE);
		uint idx = 2 * lid - (lid & (stride - 1));
		sortLocal(&local_buffer[idx], &local_buffer[idx + stride], dir);
	}

	// sync and write back
	barrier(CLK_LOCAL_MEM_FENCE);
	outArray[index].src = local_buffer[lid].src;
	outArray[index].dest = local_buffer[lid].dest;
	outArray[index].weight= local_buffer[lid].weight;
	outArray[index + MAX_LOCAL_SIZE].src = local_buffer[lid + MAX_LOCAL_SIZE].src;
	outArray[index + MAX_LOCAL_SIZE].dest = local_buffer[lid + MAX_LOCAL_SIZE].dest;
	outArray[index + MAX_LOCAL_SIZE].weight = local_buffer[lid + MAX_LOCAL_SIZE].weight;
}

__kernel void Sort_BitonicMergesortLocal(__global node* data, const uint size, const uint blocksize, uint stride)
{
	// This Kernel is basically the same as Sort_BitonicMergesortStart except of the "unrolled" part and the provided parameters
	__local node local_buffer[2 * MAX_LOCAL_SIZE];
	uint gid = get_global_id(0);
	uint groupId = get_group_id(0);
	uint lid = get_local_id(0);
	uint clampedGID = gid & (size / 2 - 1);

	uint index = groupId * (MAX_LOCAL_SIZE * 2) + lid;
	//load into local mem
	local_buffer[lid].src = data[index].src;
	local_buffer[lid].dest = data[index].dest;
	local_buffer[lid].weight= data[index].weight;
	local_buffer[lid + MAX_LOCAL_SIZE].src = data[index + MAX_LOCAL_SIZE].src;
		local_buffer[lid + MAX_LOCAL_SIZE].dest = data[index + MAX_LOCAL_SIZE].dest;
			local_buffer[lid + MAX_LOCAL_SIZE].weight = data[index + MAX_LOCAL_SIZE].weight;
	// bitonic merge
	char dir = (clampedGID & (blocksize / 2)) == 0; //same as above, % calc
#pragma unroll
	for (; stride > 0; stride >>= 1) {
		barrier(CLK_LOCAL_MEM_FENCE);
		uint idx = 2 * lid - (lid & (stride - 1));
		sortLocal(&local_buffer[idx], &local_buffer[idx + stride], dir);
	}

	// sync and write back
	barrier(CLK_LOCAL_MEM_FENCE);
	data[index].src = local_buffer[lid].src;
	data[index].dest = local_buffer[lid].dest;
	data[index].weight = local_buffer[lid].weight;
	data[index + MAX_LOCAL_SIZE].src = local_buffer[lid + MAX_LOCAL_SIZE].src;
	data[index + MAX_LOCAL_SIZE].dest = local_buffer[lid + MAX_LOCAL_SIZE].dest;
	data[index + MAX_LOCAL_SIZE].weight = local_buffer[lid + MAX_LOCAL_SIZE].weight;
}

__kernel void Sort_BitonicMergesortGlobal(__global node* data, const uint size, const uint blocksize, const uint stride)
{
	// TO DO: Kernel implementation
	uint gid = get_global_id(0);
	uint clampedGID = gid & (size / 2 - 1);

	//calculate index and dir like above
	uint index = 2 * clampedGID - (clampedGID & (stride - 1));
	char dir = (clampedGID & (blocksize / 2)) == 0; //same as above, % calc

	//bitonic merge
	node left;
	left.src= data[index].src;
	left.dest = data[index].dest;
	left.weight = data[index].weight;
	node right ;
	right.src= data[index + stride].src;
		right.dest= data[index + stride].dest;
			right.weight= data[index + stride].weight;

	sort(&left, &right, dir);

	// writeback
	data[index].src = left.src;
	data[index].dest = left.dest;
	data[index].weight = left.weight;
	data[index + stride].src = right.src;
data[index + stride].dest = right.dest;
data[index + stride].weight = right.weight;

}

int root(int v,int p[]){
 
    while(p[v] != v)
        {v = p[v];}
         
return v;
}
 
void union_ij(int i,int j,int p[]){
    if(j > i)
        p[j] = i;
    else
        p[i] = j;
}

__kernel void kruskal_algo(__global node* inArray, __global node* outArray,__global uint * p,__global uint* indicate,__global uint* test,__global uint* flag)
{	int size=2;
	

	int gid = get_global_id(0);
	
	__local uint local_flag[1];
	local_flag[0]=flag[0];
	if(size-gid-1>=0 && (size-gid-1) <size)
	{int v;
	int i;
	int k=0;
	while(atomic_cmpxchg(&flag[size-gid],1,1)!=1 && atomic_cmpxchg(&indicate[size-gid-1],indicate[size-gid-1],indicate[size-gid-1])!=2 && k<10)
	{k++;
	//barrier(CLK_LOCAL_MEM_FENCE);
	//test[size-gid-1]=atomic_cmpxchg(flag,gid,gid);

				 v=inArray[size-gid-1].src;	
				
		 
	}
	//barrier(CLK_LOCAL_MEM_FENCE);
	test[size-gid-1]=atomic_cmpxchg(&flag[size-gid-1],1,1);
	if(atomic_cmpxchg(&flag[size-gid],1,1)==1)
	{
		if(inArray[size-gid-1].weight != 999 && atomic_cmpxchg(&indicate[size-gid-1],indicate[size-gid-1],indicate[size-gid-1])!=2)
			{ v=inArray[size-gid-1].src;	
		
			while(p[v] != v)
					{v = p[v];}

			   i = v;
			   v=inArray[size-gid-1].dest;	
			while(p[v] != v)
					{v = p[v];}
			   int  j = v;
				if (i != j)
				{	atomic_xchg(&indicate[size-gid-1],1);
				indicate[size-gid-1]=1;
					outArray[size-gid-1]=inArray[size-gid-1];
					if(j > i)
						 p[j] = i;
					else
						p[i] = j;
					//test[1]=20;
				}
				else
				 {atomic_xchg(&indicate[size-gid-1],2);	
				 indicate[size-gid-1]=2;
					outArray[size-gid-1]=inArray[size-gid-1];
				//test[1]=inArray[size-gid-1].src;
				//test[0]=inArray[size-gid-1].dest;
				
				}
			
			}
	
	}
	atomic_xchg(&flag[size-gid-1],1);
	//atomic_xchg(local_flag,local_flag[0]+1);
		//test[size-gid-1]=atomic_cmpxchg(local_flag,local_flag[0],local_flag[0]);
	mem_fence(CLK_GLOBAL_MEM_FENCE);
	
	//barrier(CLK_LOCAL_MEM_FENCE);
	
	}
	
}

__kernel void test_kruskal(__global node* inArray,__global uint * p,__global uint* indicate,__global uint* test,__global uint* size,__global node* outArray)
{	//int size=2;

	int gid = get_global_id(0);
	int v;
	int i;
	if( indicate[gid]!=2){
	//if(gid>size[0] && indicate[gid]!=2){
	v=inArray[gid].src;	
				while(p[v] != v)
						{v = p[v];}

				    i = v;
				   v=inArray[gid].dest;	
				while(p[v] != v)
						{v = p[v];}
				   int  j = v;
				   if(gid==5){
				  atomic_xchg(&test[0],i);
				  atomic_xchg(&test[1],j);
	}
						outArray[gid].src=i;
						outArray[gid].dest=j;
				   if(i==j)
				   {
						indicate[gid]=2;

		   
				   }
}
}