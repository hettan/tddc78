#include "particle_list.hh"


ParticleList::ParticleList(){
  first = nullptr;
  arr = nullptr;
  size = 0;
}

//Free all particles on destroy
ParticleList::~ParticleList(){
  particle_t* next;

  while(first){
    next = first->next;
    free(first);
    first = next;
  }

  if(arr != nullptr){
    free(arr);
  }
}

void ParticleList::create_random_particles(const long n,
					   const int box_width,
					   const int box_height)
{
  //Add process_id to random seed for better random between the processes
  const int myid = MPI::COMM_WORLD.Get_rank();  
  srand(time(NULL)+myid);
  
  for(int i=0; i<n; i++){
    particle_t* p = static_cast<particle_t*>(malloc(sizeof(particle_t)));
    
    p->pcord.x = f_rand(0, box_width);
    p->pcord.y = f_rand(0, box_height);
    
    const double r = f_rand(0, MAX_INITIAL_VELOCITY);
    const double v = f_rand(0, 2*PI);
    p->pcord.vx = r*cos(v);
    p->pcord.vy = r*sin(v);

    insert(p);
  }
}

void ParticleList::insert(particle_t* p){
  if(first != nullptr){
    first->prev = p;
    p->next = first;
  }
  else
    p->next = nullptr;
  
  p->prev = nullptr;
  first = p;
  size++;
}

particle_t* ParticleList::remove(particle_t* p){
  if(p == nullptr)
    return nullptr;
  
  if(p == first)
    first = first->next;

  if(p->prev != nullptr)
    (p->prev)->next = p->next;
  
  if(p->next != nullptr)
    (p->next)->prev = p->prev;

  p->next = nullptr;
  p->prev = nullptr;

  size--;
  return p;
}

void ParticleList::erase(particle_t* p){
  remove(p);
  free(p);
}

particle_t* ParticleList::pop(){
  return remove(first);
}

void ParticleList::insert_all(ParticleList* l){
  particle_t* new_first = l->get_first();
  particle_t* last_p = new_first;

  //l is empty
  if(last_p == nullptr){
    return;
  }

  const int add_count = l->get_size();
  if(first != nullptr){
    //move to last particle in list l
    while(last_p->next != nullptr)
      last_p = last_p->next;
    
    //Link the lists together
    last_p->next = first;
    first->prev = last_p;
  }
  
  //update new first and empty the other list
  first = new_first;
  l->clear();
  size += add_count;
}

//Returns a sequencial memory space with all particles
particle_t* ParticleList::get_arr_copy(){
  //allocate the sequencial memory if not done already
  if(arr == nullptr)
    arr = (particle_t*)malloc(MAX_SEND_BUFF_SIZE*sizeof(particle_t));
  
  particle_t* p = first;
  const int count = size;
  for(int i=0; i<count; i++){
    arr[i] = *p;
    p = p->next;
  }

  return arr;
}

//Generate random values for n particles
double ParticleList::f_rand(const double f_min, const double f_max)
{
  const double f = (double)rand() / RAND_MAX;
  return f_min + f * (f_max - f_min);
}

int get_size(particle_t* arr)
{
  int counter=0;
  while(arr != nullptr){    
    counter++;
    arr=arr->next;
  }
  
  return counter;
}

/*
int main(){
  ParticleList* a = new ParticleList();
  a->create_random_particles(10, 100, 100);
  cout << a->get_size() << endl;
  cout << get_size(a->get_first()) << endl;

  for(int i=0; i<5; i++)
    a->pop();  
  cout << a->get_size() << endl;
  cout << get_size(a->get_first()) << endl;
  
  cout << a << endl;
  particle_t* p = a->get_first();
  for(int i=0; i<3; i++)
    p = p->next;

  a->erase(p);
  cout << a << endl;

  cout << a->get_size() << endl;
  particle_t* arr = a->get_arr_copy(); 
  const int size = a->get_size();
  for(int i=0; i < size; i++)
    cout << arr+i << endl;
  
  particle_t* tmp;
  p = a->get_first();
  while(p){
    tmp = p->next;
    a->erase(p);
    p=tmp;
  }
  
  cout << a->get_size() << endl;
  //particle_t* arr = a->get_arr_copy(); 
  //const int size = a->get_size();
  for(int i=0; i < size; i++)
    cout << arr+i << endl;
  

  return 0;
}
*/
