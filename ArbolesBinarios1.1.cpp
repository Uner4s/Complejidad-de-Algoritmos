/*
	Archivo de Complejidad.
    Autor: Nicolas Fernando Perez Poblete.
    Profesor: Oscar Rojas.
*/

#include <iostream>
#include <vector>
#include <map>
#define lonCad 15
#include <fstream>
#include <string>
#include <cstdio>
#include <sys/time.h>
#include <time.h>
#include <cstdlib>
#include <string.h>

using namespace std;

ofstream Archivo_abb_datos;
ofstream Archivo_abb_freq;
ofstream Archivo_avl_datos;
ofstream Archivo_avl_freq;
ofstream Archivo_23_datos;
ofstream Archivo_23_freq;
ofstream Archivo_b_datos;
ofstream Archivo_b_freq;
int cambio_altura; //Cambios de altura
int auxavl;
int altura = 0; //Altura desde los padres a las hijas
int global_suma=0;//sumas relevantes
int global_resta=0;//restas relevantes
int global_producto=0;//productos relevantes
int complejidad_ABB=0;//lamado a funciones
int complejidad_AVL=0;//llamado a funciones
int complejidad_B=0;//llamado a funciones

struct timespec tiempo1, tiempo2, latencia;

//*************************************************************************************************************************
void diff_time( timespec *t_fin, timespec *t_ini, timespec *delta )
{
    if( ( (*t_fin).tv_nsec - (*t_ini).tv_nsec ) < 0 )
    {
        if( (*t_fin).tv_sec == (*t_ini).tv_sec )
        {
            (*delta).tv_sec  = 0;
            (*delta).tv_nsec = 1000000000 + (*t_fin).tv_nsec - (*t_ini).tv_nsec;
        }
        else
        {
            (*delta).tv_sec  = (*t_fin).tv_sec - (*t_ini).tv_sec - 1;
            (*delta).tv_nsec = 1000000000 + (*t_fin).tv_nsec - (*t_ini).tv_nsec;
        }
    }
    else
    {
        if( (*t_fin).tv_sec == (*t_ini).tv_sec )
        {
            (*delta).tv_sec  = 0;
            (*delta).tv_nsec = (*t_fin).tv_nsec - (*t_ini).tv_nsec;
        }
        else
        {
            (*delta).tv_sec  = (*t_fin).tv_sec - (*t_ini).tv_sec;
            (*delta).tv_nsec = (*t_fin).tv_nsec - (*t_ini).tv_nsec;
        }
    }
}
//*************************************************************************************************************************
typedef struct puntero{
     int Hoja_abb;
     int frecuencia;
     struct puntero *izquierda, *derecha;
}*ABB; 
 
ABB crearNodo(int x, int y) // Crea un nuevo nodo en el arbol
{
     ABB nuevoNodo = new(struct puntero);
     nuevoNodo->Hoja_abb = x;
     nuevoNodo->izquierda = NULL;
     nuevoNodo->derecha = NULL;
     nuevoNodo->frecuencia = 0;
     complejidad_ABB++;// Se va calculando la complejidad de computo.
 
     return nuevoNodo;
}
void insertar(ABB &arbol, int x, int y) // Inserta un numero mandado al arbol
{
     if(arbol==NULL) {
           arbol = crearNodo(x,y); // crea un nuevo nodo
           complejidad_ABB++; // Se va calculando la complejidad de computo.
     }
     else if(x==arbol->Hoja_abb){
     	arbol->frecuencia=y;
	 }
     else if(x < arbol->Hoja_abb){
     	  complejidad_ABB++;// Se va calculando la complejidad de computo.
          insertar(arbol->izquierda, x, y);
      }
     else if(x > arbol->Hoja_abb){
     	  complejidad_ABB++;// Se va calculando la complejidad de computo.	
          insertar(arbol->derecha, x, y);
      }   
}

bool ABB_Tree_Busqueda(ABB arbol, int dato) // BUSQUEDA ABB
{
     int Booleano=0;   // 0 indica que lo encontre
     if(arbol==NULL)
        return Booleano;
     if(dato<arbol->Hoja_abb)
         Booleano = ABB_Tree_Busqueda(arbol->izquierda, dato);
     else if(dato> arbol->Hoja_abb)
         Booleano = ABB_Tree_Busqueda(arbol->derecha, dato);
     else
        Booleano = 1;   // son iguales, lo encontre
     return Booleano;
}
//*************************************************************************************************************************
//*************************************************************************************************************************
//*************************************************************************************************************************
//*************************************************************************************************************************
typedef struct Nodo_AVL{
    int clave;
    int bal; /* Factor de balance  */
    int frecuencia;
    struct Nodo_AVL *left, *right;
} nodo, *arbolAVL;

static arbolAVL Rotacion_AVL_Izquierda(arbolAVL Arbol_Auxiliar){
    arbolAVL temp;
    int x,y;
    temp = Arbol_Auxiliar;
    Arbol_Auxiliar = Arbol_Auxiliar->right;
    temp->right = Arbol_Auxiliar->left;
    Arbol_Auxiliar->left = temp;
    x = temp->bal; 
    y = Arbol_Auxiliar->bal; 
    temp->bal = x-1-max(y, 0);
    global_resta++;
    Arbol_Auxiliar->bal = min(x-2+min(y, 0), y-1);
    global_resta++;
    global_suma++;
    complejidad_AVL++;
    return Arbol_Auxiliar;
}

static arbolAVL Rotacion_AVL_Derecha(arbolAVL Arbol_Auxiliar){
    arbolAVL temp = Arbol_Auxiliar;
    int x,y;
    Arbol_Auxiliar = Arbol_Auxiliar->left;
    temp->left = Arbol_Auxiliar->right;
    Arbol_Auxiliar->right = temp;
    x = temp->bal;
    y = Arbol_Auxiliar->bal;   
    temp->bal = x+1-min(y, 0);  
    global_suma++;
    Arbol_Auxiliar->bal = max(x+2+max(y, 0), y+1); 
    global_suma++;
    complejidad_AVL++;
    return Arbol_Auxiliar;
}

int Altura(void){
	complejidad_AVL++;
    return altura;
}

arbolAVL CreaNodo(int auxavl, int frec){
    arbolAVL Arbol_Auxiliar;
    Arbol_Auxiliar = (arbolAVL) malloc(sizeof(nodo));
    Arbol_Auxiliar->clave=auxavl;
    Arbol_Auxiliar->left=0;
    Arbol_Auxiliar->right=0;
    Arbol_Auxiliar->frecuencia = 0;
    complejidad_AVL++;
    return Arbol_Auxiliar;
}

arbolAVL Insertar_y_Balancear(arbolAVL Arbol_Auxiliar, int frec)
{
    if (Arbol_Auxiliar == NULL){
        Arbol_Auxiliar = CreaNodo(auxavl, frec);
        Arbol_Auxiliar->bal = 0; /* Los dos hijos son nulos */
        cambio_altura = 1; /* Marca necesidad de revisar balances */
        complejidad_AVL++;// Se va calculando la complejidad de computo.
        return Arbol_Auxiliar; /* retorna puntero al insertado */
    }
    else if (Arbol_Auxiliar->clave < auxavl){
    	complejidad_AVL++;// Se va calculando la complejidad de computo.
        //desciende por la derecha
        Arbol_Auxiliar->right = Insertar_y_Balancear(Arbol_Auxiliar->right, frec);
        //se pasa por la siguiente linea en la revision ascendente
        Arbol_Auxiliar->bal += cambio_altura; /* Incrementa factor de balance */
        global_suma++;
    }
    else if (Arbol_Auxiliar->clave > auxavl){
    	complejidad_AVL++;// Se va calculando la complejidad de computo.
        //desciende por la izquierda
        Arbol_Auxiliar->left = Insertar_y_Balancear(Arbol_Auxiliar->left, frec);
        //se corrige en el ascenso
        Arbol_Auxiliar->bal -= cambio_altura; /* Decrementa balance */
        global_resta++;
    }
    else  
    {
    	Arbol_Auxiliar->frecuencia = frec;
    	complejidad_AVL++;// Se va calculando la complejidad de computo.
        cambio_altura = 0;
    }
    if (cambio_altura == 0) /* No hay que rebalancear. Sigue el ascenso */
        return Arbol_Auxiliar;
    if (Arbol_Auxiliar->bal < -1){
        /* esta desbalanceado por la izquierda.*/
        if (Arbol_Auxiliar->left->bal > 0){
        	complejidad_AVL++;// Se va calculando la complejidad de computo.
            Arbol_Auxiliar->left = Rotacion_AVL_Izquierda(Arbol_Auxiliar->left);
        }
        complejidad_AVL++;// Se va calculando la complejidad de computo.
        Arbol_Auxiliar = Rotacion_AVL_Derecha(Arbol_Auxiliar);
        cambio_altura = 0; /* El subarbol no aumenta su altura */
    }
    else if (Arbol_Auxiliar->bal > 1){
        /* Si queda desbalanceado por la derecha.*/
        if (Arbol_Auxiliar->right->bal < 0){
        	complejidad_AVL++;// Se va calculando la complejidad de computo.
            /* Si hijo derecho esta cargado a la izquierda Caso d.*/
            Arbol_Auxiliar->right = Rotacion_AVL_Derecha(Arbol_Auxiliar->right);
        }
        Arbol_Auxiliar = Rotacion_AVL_Izquierda(Arbol_Auxiliar); /* Caso c.*/
        cambio_altura = 0; /* El subarbol no aumenta su altura */
        complejidad_AVL++;// Se va calculando la complejidad de computo.
    }
    else if (Arbol_Auxiliar->bal == 0){/* La inserccion lo balanceo */
    	complejidad_AVL++;// Se va calculando la complejidad de computo.
        cambio_altura = 0; /* El subarbol no aumenta su altura. Caso a*/
    }
    else {
        complejidad_AVL++;// Se va calculando la complejidad de computo.
		cambio_altura = 1; /* Propaga ascendentemente la necesidad de rebalancear */
    }
    complejidad_AVL++;// Se va calculando la complejidad de computo.
    return Arbol_Auxiliar;
}

arbolAVL InsertarAVL(int clave, arbolAVL Avl, int frec){
	complejidad_AVL++;
    auxavl = clave; 
    Avl = Insertar_y_Balancear(Avl, frec);
    if (cambio_altura == 1){
    	altura++;
    	global_suma++;
	}
    return Avl;
}

void Inorden_AVL_Arbol(arbolAVL Arbol_Auxiliar){
    if (Arbol_Auxiliar != NULL){
        Inorden_AVL_Arbol(Arbol_Auxiliar->left);
        Archivo_avl_datos<<Arbol_Auxiliar->clave<<" ";
		Archivo_avl_freq<<Arbol_Auxiliar->frecuencia<<" ";
        Inorden_AVL_Arbol(Arbol_Auxiliar->right);
    }
}

int AVL_Tree_Busqueda(arbolAVL Arbol_Auxiliar, int dato){ // BUSQUEDA AVL
    int Booleano=0;   // 0 indica que lo encontre
    if(Arbol_Auxiliar==NULL)
    	return Booleano;
    if(dato<Arbol_Auxiliar->clave)
    	Booleano = AVL_Tree_Busqueda(Arbol_Auxiliar->left, dato);
    else if(dato> Arbol_Auxiliar->clave)
    	Booleano = AVL_Tree_Busqueda(Arbol_Auxiliar->right, dato);
    else if(dato == Arbol_Auxiliar->clave)
    	Booleano = 1;  
    return Booleano;
}

void Inorden_ABB_Arbol(ABB arbol){ //Funcion que guarda en archivo los datos de los arboles abb.
    if(arbol!=NULL){
        Inorden_ABB_Arbol(arbol->izquierda);
        Archivo_abb_datos<<arbol->Hoja_abb << " "; 
		Archivo_abb_freq<<arbol->frecuencia<<" ";
        Inorden_ABB_Arbol(arbol->derecha);
    }
}
//**********************************************************************************************************************************
class Nodo_Arbol_B{
    int *Hojas_B;  // An array of Hojas_B
    int *frecuencias; // Arreglo de frecuencias
    int t;      // Minimum degree (defines the range for number of Hojas_B)
    Nodo_Arbol_B **C; // An array of child pointers
    int n;     // Current number of Hojas_B
    bool leaf; // Is true when node is leaf. Otherwise false
 
public:
    Nodo_Arbol_B(int _t, bool _leaf);   // Constructor
	void Recorrer_Arbol_23();
    void Recorrer_Arbol_B();
    Nodo_Arbol_B *search(int k);
	Nodo_Arbol_B *AgregarFreq(int k, int freq); 
    int findKey(int k);
    void insertNonFull(int k);
    void splitChild(int i, Nodo_Arbol_B *y);
    int getPred(int idx);
    int getSucc(int idx);
    void fill(int idx);
    void borrowFromPrev(int idx);
    void borrowFromNext(int idx);
    void merge(int idx);
    friend class Arbol_B;
};
 
class Arbol_B{
    int t;
public:
    Nodo_Arbol_B *root;
    Arbol_B(int _t){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        root = NULL;
        t = _t;
    }
    void Recorrer_Arbol_23(){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        if (root != NULL) root->Recorrer_Arbol_23();
    }
    void Recorrer_Arbol_B(){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        if (root != NULL) root->Recorrer_Arbol_B();
    }
    Nodo_Arbol_B* search(int k){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        return (root == NULL)? NULL : root->search(k);
    }
    Nodo_Arbol_B* AgregarFreq(int k, int freq){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        return (root == NULL)? NULL : root->AgregarFreq(k,freq);
    }  
    void insert(int k);
 	void VolverNULLpuntero();
    //void remove(int k);
 
};
typedef class Nodo_Arbol_B * B_23;
 
Nodo_Arbol_B::Nodo_Arbol_B(int t1, bool leaf1){
    // Copy the given minimum degree and leaf property
    t = t1;
    leaf = leaf1;
 	complejidad_B++;// Se va calculando la complejidad de computo.
    // Allocate memory for maximum number of possible Hojas_B
    // and child pointers
    Hojas_B = new int[2*t-1];
    global_producto=global_producto+3;
    global_resta=global_resta+3;
    frecuencias = new int[2*t-1];
    C = new Nodo_Arbol_B *[2*t];
    n = 0;
}

int Nodo_Arbol_B::findKey(int k){
    int idx=0;
    complejidad_B++;
    while (idx<n && Hojas_B[idx] < k){
    	global_suma++;// Se va calculando la complejidad de computo.
        ++idx;
    }
    return idx;
}
 
int Nodo_Arbol_B::getPred(int idx){
    // Keep moving to the right most node until we reach a leaf
    Nodo_Arbol_B *cur=C[idx];
    while (!cur->leaf){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        cur = cur->C[cur->n];
    }
    // Return the last key of the leaf
    global_resta++;
    return cur->Hojas_B[cur->n-1];
}
 
int Nodo_Arbol_B::getSucc(int idx){
    // Keep moving the left most node starting from C[idx+1] until we reach a leaf
    Nodo_Arbol_B *cur = C[idx+1];
    global_suma++;
    complejidad_B++;// Se va calculando la complejidad de computo.
    while (!cur->leaf){
        cur = cur->C[0];
    }
    // Return the first key of the leaf
    complejidad_B++;// Se va calculando la complejidad de computo.
    return cur->Hojas_B[0];
}
 
void Nodo_Arbol_B::fill(int idx){
 	complejidad_B++;// Se va calculando la complejidad de computo.
    if (idx!=0 && C[idx-1]->n>=t){
    	global_resta++;
    	complejidad_B++;// Se va calculando la complejidad de computo.
        borrowFromPrev(idx);
    }
    // If the next child(C[idx+1]) has more than t-1 Hojas_B, borrow a key
    // from that child
    else if (idx!=n && C[idx+1]->n>=t){
    	global_suma++;
    	complejidad_B++;// Se va calculando la complejidad de computo.
        borrowFromNext(idx);
    }
    else{
        if (idx != n){
        	complejidad_B++;// Se va calculando la complejidad de computo.
            merge(idx);
        }
        else{
        	complejidad_B++;// Se va calculando la complejidad de computo.
            merge(idx-1);
        }
    }
    return;
}
 
void Nodo_Arbol_B::borrowFromPrev(int idx){
 	complejidad_B++;// Se va calculando la complejidad de computo.
    Nodo_Arbol_B *child=C[idx];
    Nodo_Arbol_B *sibling=C[idx-1];
    global_resta++;
    for (int i=child->n-1; i>=0; --i){
    	global_resta++;
    	complejidad_B++;// Se va calculando la complejidad de computo.
        child->Hojas_B[i+1] = child->Hojas_B[i];
        global_suma++;
    }
    // If C[idx] is not a leaf, move all its child pointers one step ahead
    if (!child->leaf){	
    	complejidad_B++;
        for(int i=child->n; i>=0; --i){
        	complejidad_B++;// Se va calculando la complejidad de computo.
            child->C[i+1] = child->C[i];
            global_suma++;
        }
    }
    // Setting child's first key equal to Hojas_B[idx-1] from the current node
    child->Hojas_B[0] = Hojas_B[idx-1];
    global_resta++;
 	complejidad_B++;// Se va calculando la complejidad de computo.
    // Moving sibling's last child as C[idx]'s first child
    if (!leaf){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        child->C[0] = sibling->C[sibling->n];
    }
    complejidad_B++;
    Hojas_B[idx-1] = sibling->Hojas_B[sibling->n-1];
    child->n += 1;
    sibling->n -= 1;
    global_resta++;
    global_suma++;
    return;
}

void Nodo_Arbol_B::borrowFromNext(int idx){
    Nodo_Arbol_B *child=C[idx];
    Nodo_Arbol_B *sibling=C[idx+1];
 	complejidad_B++;// Se va calculando la complejidad de computo.
    // Hojas_B[idx] is inserted as the last key in C[idx]
    child->Hojas_B[(child->n)] = Hojas_B[idx];
    // Sibling's first child is inserted as the last child
    // into C[idx]
    if (!(child->leaf)){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        child->C[(child->n)+1] = sibling->C[0];
        global_suma++;
    }
    //The first key from sibling is inserted into Hojas_B[idx]
    Hojas_B[idx] = sibling->Hojas_B[0];
    // Moving all Hojas_B in sibling one step behind
    for (int i=1; i<sibling->n; ++i){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        sibling->Hojas_B[i-1] = sibling->Hojas_B[i];
    }
    // Moving the child pointers one step behind
    if (!sibling->leaf){
        for(int i=1; i<=sibling->n; ++i){
        	complejidad_B++;// Se va calculando la complejidad de computo.
            sibling->C[i-1] = sibling->C[i];
        }
    }
    // Increasing and decreasing the key count of C[idx] and C[idx+1]
    // respectively
    child->n += 1;
    global_suma++;
    sibling->n -= 1;
 	complejidad_B++;// Se va calculando la complejidad de computo.
    return;
}

void Nodo_Arbol_B::merge(int idx){
    Nodo_Arbol_B *child = C[idx];
    Nodo_Arbol_B *sibling = C[idx+1];
    // Pulling a key from the current node and inserting it into (t-1)th
    // position of C[idx]
    child->Hojas_B[t-1] = Hojas_B[idx];
    // Copying the Hojas_B from C[idx+1] to C[idx] at the end
    for (int i=0; i<sibling->n; ++i){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        child->Hojas_B[i+t] = sibling->Hojas_B[i];
    }
    // Copying the child pointers from C[idx+1] to C[idx]
    if (!child->leaf){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        for(int i=0; i<=sibling->n; ++i){
        	complejidad_B++;// Se va calculando la complejidad de computo.
            child->C[i+t] = sibling->C[i];
        }
    }
    // Moving all Hojas_B after idx in the current node one step before -
    // to fill the gap created by moving Hojas_B[idx] to C[idx]
    for (int i=idx+1; i<n; ++i){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        Hojas_B[i-1] = Hojas_B[i];
    }
    // Moving the child pointers after (idx+1) in the current node one
    // step before
    for (int i=idx+2; i<=n; ++i){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        C[i-1] = C[i];
    }
    // Updating the key count of child and the current node
    child->n += sibling->n+1;
    n--;
 	complejidad_B++;// Se va calculando la complejidad de computo.
    // Freeing the memory occupied by sibling
    delete(sibling);
    return;
}
 
void Arbol_B::insert(int k)
{
    // If tree is empty
    if (root == NULL)
    {
    	complejidad_B++;// Se va calculando la complejidad de computo.
        // Allocate memory for root
        root = new Nodo_Arbol_B(t, true);
        root->Hojas_B[0] = k;  // Insert key
        root->n = 1;  // Update number of Hojas_B in root
    }
    else // If tree is not empty
    {
        // If root is full, then tree grows in height
        if (root->n == 2*t-1){
        	complejidad_B++;// Se va calculando la complejidad de computo.
            // Allocate memory for new root
            Nodo_Arbol_B *s = new Nodo_Arbol_B(t, false);
            // Make old root as child of new root
            s->C[0] = root;
            // Split the old root and move 1 key to the new root
            s->splitChild(0, root);
            // New root has two children now.  Decide which of the
            // two children is going to have new key
            int i = 0;
            if (s->Hojas_B[0] < k)
                i++;
            global_suma++;    
            s->C[i]->insertNonFull(k);
            complejidad_B++;// Se va calculando la complejidad de computo.
            // Change root
            root = s;
        }
        else{  // If root is not full, call insertNonFull for root
        	complejidad_B++;// Se va calculando la complejidad de computo.
            root->insertNonFull(k);
        }
    }
}
 
void Nodo_Arbol_B::insertNonFull(int k){
    // Initialize index as index of rightmost element
    int i = n-1;
    // If this is a leaf node
    if (leaf == true)
    {
    	complejidad_B++;// Se va calculando la complejidad de computo.
        // The following loop does two things
        // a) Finds the location of new key to be inserted
        // b) Moves all greater Hojas_B to one place ahead
        complejidad_B++;// Se va calculando la complejidad de computo.
        while (i >= 0 && Hojas_B[i] > k)
        {
        	global_suma++;
            Hojas_B[i+1] = Hojas_B[i];
            i--;
            global_resta++;
        }
 
        // Insert the new key at found location
        Hojas_B[i+1] = k;
        n = n+1;
    }
    else // If this node is not leaf
    {
        // Find the child which is going to have the new key
        while (i >= 0 && Hojas_B[i] > k){
        	complejidad_B++;// Se va calculando la complejidad de computo.
            i--;
        }
        // See if the found child is full
        if (C[i+1]->n == 2*t-1){
        	complejidad_B++;// Se va calculando la complejidad de computo.
            // If the child is full, then split it
            splitChild(i+1, C[i+1]);
 			global_suma++;
			global_suma++;      
            // After split, the middle key of C[i] goes up and
            // C[i] is splitted into two.  See which of the two
            // is going to have the new key
            if (Hojas_B[i+1] < k)
                i++;
        }
        complejidad_B++;// Se va calculando la complejidad de computo.
        C[i+1]->insertNonFull(k);
        global_suma++;   
    }
}
 
void Nodo_Arbol_B::splitChild(int i, Nodo_Arbol_B *y){
    // Create a new node which is going to store (t-1) Hojas_B
    // of y
    Nodo_Arbol_B *z = new Nodo_Arbol_B(y->t, y->leaf);
    z->n = t - 1;
 	complejidad_B++;// Se va calculando la complejidad de computo.
    // Copy the last (t-1) Hojas_B of y to z
    for (int j = 0; j < t-1; j++){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        z->Hojas_B[j] = y->Hojas_B[j+t];
        global_suma++;
    }
    // Copy the last t children of y to z
    if (y->leaf == false){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        for (int j = 0; j < t; j++){
        	complejidad_B++;// Se va calculando la complejidad de computo.
            z->C[j] = y->C[j+t];
        }
    }
    // Reduce the number of Hojas_B in y
    y->n = t - 1;
    global_resta++;
 	complejidad_B++;// Se va calculando la complejidad de computo.
    // Since this node is going to have a new child,
    // create space of new child
    for (int j = n; j >= i+1; j--){
    	complejidad_B++;// Se va calculando la complejidad de computo.
        C[j+1] = C[j];
        global_suma++;
    }
    // Link the new child to this node
    C[i+1] = z;
    // A key of y will move to this node. Find location of
    // new key and move all greater Hojas_B one space ahead
    complejidad_B++;// Se va calculando la complejidad de computo.
    for (int j = n-1; j >= i; j--){
		global_suma++;
        Hojas_B[j+1] = Hojas_B[j];
    }
    // Copy the middle key of y to this node
    Hojas_B[i] = y->Hojas_B[t-1];
    // Increment count of Hojas_B in this node
    n = n + 1;
}
 
// Function to traverse all nodes in a subtree rooted with this node
void Nodo_Arbol_B::Recorrer_Arbol_23(){
    // There are n Hojas_B and n+1 children, travers through n Hojas_B
    // and first n children
    int i;
    for (i = 0; i < n; i++){
        // If this is not leaf, then before printing key[i],
        // traverse the subtree rooted with child C[i].
        if (leaf == false)
            C[i]->Recorrer_Arbol_23();
        //cout << " " << Hojas_B[i];
        Archivo_23_datos<<Hojas_B[i]<<" ";
        Archivo_23_freq<<frecuencias[i]<<" ";
    }
    // Print the subtree rooted with last child
    if (leaf == false)
        C[i]->Recorrer_Arbol_23();
}


void Nodo_Arbol_B::Recorrer_Arbol_B(){
    // There are n Hojas_B and n+1 children, travers through n Hojas_B
    // and first n children
    int i;
    for (i = 0; i < n; i++){
        // If this is not leaf, then before printing key[i],
        // traverse the subtree rooted with child C[i].
        if (leaf == false)
            C[i]->Recorrer_Arbol_B();
        //cout << " " << Hojas_B[i];
        Archivo_b_datos<<Hojas_B[i]<<" ";
        Archivo_b_freq<<frecuencias[i]<<" ";
    }
    // Print the subtree rooted with last child
    if (leaf == false)
        C[i]->Recorrer_Arbol_B();
}


Nodo_Arbol_B *Nodo_Arbol_B::AgregarFreq(int k, int freq){
    // Find the first key greater than or equal to k
    int i = 0;
    while (i < n && k > Hojas_B[i])
        i++;
    // If the found key is equal to k, return this node
    if (Hojas_B[i] == k){
    	frecuencias[i]=freq;
    	if (Hojas_B[i] == k)
        	return this;
    }
    // If key is not found here and this is a leaf node
    if (leaf == true)
        return NULL;
    // Go to the appropriate child
    return C[i]->AgregarFreq(k,freq);
}

// Function to search key k in subtree rooted with this node
Nodo_Arbol_B *Nodo_Arbol_B::search(int k){
    // Find the first key greater than or equal to k
    int i = 0;
    while (i < n && k > Hojas_B[i])
        i++;
    // If the found key is equal to k, return this node
    if (Hojas_B[i] == k){
    	//cout<<"si esta"<<endl;
        return this;
    }
    // If key is not found here and this is a leaf node
    if (leaf == true)
        return NULL;
    // Go to the appropriate child
    return C[i]->search(k);
}

void Arbol_B::VolverNULLpuntero(){
	root = NULL;
	return;
}
//*************************************************************************************************************************
//*************************************************************************************************************************
int esta(int numero, vector<int> UV, int i2) //Funcion que sirve para crear un vector con los elementos una sola vez.
{
	int cont=0;
	int respuesta=0;
	while(cont<i2){
		if(numero==UV[cont]){
			respuesta++;
			cont++;
		}
		else{
			cont++;
		}
	}
	if(respuesta>0)
		return 1;
	else
		return 0;
}

int frecuencia(int numero, vector<int> elementos, int cantidad)// Esta funcion retorna las frecuencias de un numero especifico
{
	int i=0;
	int freq=0;	
	while(i<cantidad){
		if(numero==elementos[i]){
			freq++;
			i++;
		}
		else{
			i++;
		}
	}
	return freq;	
}
//*************************************************************************************************************************
//*************************************************************************************************************************
void Crear_Arboles(map <int , vector<int> > Mapa_Datos, map <int , vector<int> > Mapa_Frecuencias, map <int , vector<int> >	Busqueda, map <int , vector<int> >	Map_Datos ){
	ofstream complejidad_temporal;
	complejidad_temporal.open("ComplejidadT.dat");
    cout<<"Comenzando a crear los arboles"<<endl;
    cout<<endl;
	map <int, vector < int> >:: iterator iterador;
	map <int, vector < int> >:: iterator iteradorfreq;
    map <int , ABB> Mapa_ABB_Datos; // Mapa de arbol abb de los datos
    map <int , arbolAVL> Mapa_AVL_Datos; // Mapa de arbol avl de los datos
    map <int , B_23> Mapa_23_Datos; // Mapa de arbol 2-3 de los datos
    map <int , B_23> Mapa_BTREE_Datos; // Mapa de arbol b de los datos
	map <int , vector<int> > Complejidad_abb;
	map <int , vector<int> > Complejidad_avl;
	map <int , vector<int> > Complejidad_23;
	map <int , vector<int> > Complejidad_b;
	vector<int> operacion_abb;
    vector<int> operacion_avl;
    vector<int> operacion_23;
    vector<int> operacion_b;
    ABB ABB_Aux;
    arbolAVL AVL_Aux; // Variables auxiliares.
    B_23 B_23_Aux;
	ABB Arbol_ABB_D;
    ABB Arbol_ABB_F;
	Arbol_ABB_D = NULL;
    Arbol_ABB_F = NULL;
    arbolAVL Arbol_AVL_D;
    arbolAVL Arbol_AVL_F;
    Arbol_AVL_D = NULL;
    Arbol_AVL_F = NULL;
    long int tiempoComputo;
	Archivo_abb_datos.open("Arbol_ABB_Datos.dat");
	Archivo_abb_freq.open("Arbol_ABB_Frecuencias.dat");
//------------------------------------------------------------------------------------------------------------
	int key;
	complejidad_temporal<<"Arbol ABB"<<endl;
	iteradorfreq = Mapa_Datos.begin();
    for(iterador =Map_Datos.begin(); iterador != Map_Datos.end(); ++iterador){ // Se recorre con un iterador el mapa de datos del archivo.
    	key=iterador->first;
    	Archivo_abb_datos<<"Arbol ABB datos ()KEY: "<<key<<endl;
    	Archivo_abb_freq<<"Arbol ABB frecuencias ()KEY: "<<key<<endl;
    	clock_gettime(CLOCK_MONOTONIC, &tiempo1);
        for (int i =0;i<iterador->second.size(); i++){ // hasta que llegue al final de la primera linea.
            insertar( Arbol_ABB_D, iterador->second[i], 0); // Va mandando datos de la fila para crear el arbol abb.
        }
        clock_gettime(CLOCK_MONOTONIC, &tiempo2);
   		diff_time(&tiempo2, &tiempo1, &latencia);
   		tiempoComputo = latencia.tv_sec*1000000000 + latencia.tv_nsec;   
		complejidad_temporal<<tiempoComputo<<" ";
        operacion_abb.push_back(complejidad_ABB);
        operacion_abb.push_back(global_suma);
        operacion_abb.push_back(global_resta);
        operacion_abb.push_back(global_producto);
        Complejidad_abb.insert(pair<int,vector<int> >(key,operacion_abb));
        complejidad_ABB=0;
        global_suma=0;
        global_resta=0;
        global_producto=0;
        operacion_abb.clear();
        for (int i =0;i<iteradorfreq->second.size(); i++){ // hasta que llegue al final de la primera linea.
            insertar( Arbol_ABB_D, iteradorfreq->second[i], Mapa_Frecuencias[key][i]); // Va mandando datos de la fila para crear el arbol abb.
        }
        ++iteradorfreq;
        Inorden_ABB_Arbol(Arbol_ABB_D); // Se llama la funcion para escribir el arbol abb en un archivo.
        Archivo_abb_datos<<endl;
        Archivo_abb_datos<<endl;
        Archivo_abb_freq<<endl;
        Archivo_abb_freq<<endl;
        Mapa_ABB_Datos.insert(pair<int, ABB  >(iterador->first,Arbol_ABB_D)); // Se ingresa el arbol creado a un mapa de arboles abb
        Arbol_ABB_D = NULL;
    }
	cout<<"Arbol ABB Creado."<<endl;    
    Archivo_abb_datos.close();
	Archivo_abb_freq.close(); 
	complejidad_temporal<<endl;
//------------------------------------------------------------------------------------------------------------
	complejidad_temporal<<"Arbol AVL"<<endl;
	Archivo_avl_datos.open("Arbol_AVL_Datos.dat"); // Se abre un archivo para guardar los arboles avl de datos
	Archivo_avl_freq.open("Arbol_AVL_Frecuencias.dat"); // Se abre un archivo para guardar los arboles avl de frecuencias
	iteradorfreq = Mapa_Datos.begin();
    for(iterador =Map_Datos.begin(); iterador != Map_Datos.end(); ++iterador){// Se recorre con un iterador el mapa de datos del archivo.
    	key=iterador->first;
    	Archivo_avl_datos<<"Arbol AVL datos ()KEY: "<<key<<endl;
    	Archivo_avl_freq<<"Arbol AVL frecuencias ()KEY: "<<key<<endl;
    	clock_gettime(CLOCK_MONOTONIC, &tiempo1);
        for (int i =0;i<iterador->second.size(); i++){// hasta que llegue al final de la primera linea.
            Arbol_AVL_D=InsertarAVL(iterador->second[i], Arbol_AVL_D, 0); // va mandando datos de la fila para crear un arbol AVL.
        }
        clock_gettime(CLOCK_MONOTONIC, &tiempo2);
   		diff_time(&tiempo2, &tiempo1, &latencia);
   		tiempoComputo = latencia.tv_sec*1000000000 + latencia.tv_nsec;
   		complejidad_temporal<<tiempoComputo<<" ";
        operacion_avl.push_back(complejidad_AVL);
        operacion_avl.push_back(global_suma);
        operacion_avl.push_back(global_resta);
        operacion_avl.push_back(global_producto);
        Complejidad_avl.insert(pair<int,vector<int> >(key,operacion_avl));
        operacion_avl.clear();
        for (int i =0;i<iteradorfreq->second.size(); i++){// hasta que llegue al final de la primera linea.
            Arbol_AVL_D=InsertarAVL(iteradorfreq->second[i], Arbol_AVL_D, Mapa_Frecuencias[key][i]); // va mandando datos de la fila para crear un arbol AVL.
        }
        ++iteradorfreq;
        complejidad_AVL=0;
        global_suma=0;
        global_resta=0;
        global_producto=0;
        Inorden_AVL_Arbol(Arbol_AVL_D); // Se llama la funcion para guardar en un archivo el arbol avl creado.
        Archivo_avl_datos<<endl;
        Archivo_avl_datos<<endl;
        Archivo_avl_freq<<endl;
        Archivo_avl_freq<<endl;
        Mapa_AVL_Datos.insert(pair<int, arbolAVL  >(iterador->first,Arbol_AVL_D)); // S van guardando los arboles avl creados en un map.
        Arbol_AVL_D = NULL; // el arbol se hace nulo.
    }
	cout<<"Arbol AVL Creado."<<endl;    
	Archivo_avl_datos.close(); // Se cierrra el archivo donde se guardan los arboles avl de datos.
    Archivo_avl_freq.close();// Se cierra el archivo donde se guardan los arboles avl de frecuencia.
    complejidad_temporal<<endl;
//------------------------------------------------------------------------------------------------------------
	complejidad_temporal<<"Arbol 2-3"<<endl;
	Archivo_23_datos.open("Arbol_23_datos.dat"); // Se abre un archivo de escritura para guardar arboles 2-3 de datos
	Archivo_23_freq.open("Arbol_23_Frecuencias.dat"); // Se abre un archivo de escritura para guardar arboles 2-3 de frecuencias.
    Arbol_B t(3);// Se hace b = 3 para simular el arbol 2-3.
    for(iterador =Mapa_Datos.begin(); iterador != Mapa_Datos.end(); ++iterador){// Se recorre con un iterador el mapa de datos del archivo.
    	key=iterador->first;
    	Archivo_23_datos<<"Arbol 2-3 datos ()KEY: "<<key<<endl;
    	Archivo_23_freq<<"Arbol 2-3 frecuencias ()KEY: "<<key<<endl;
    	
    	clock_gettime(CLOCK_MONOTONIC, &tiempo1);
        for (int i =0;i<iterador->second.size(); i++){ // Se recorre hasta el final de los datos de las key.
            t.insert(iterador->second[i]);// Se van insertando los datos de la fila en el arbol 2-3 de datos.
        }
        clock_gettime(CLOCK_MONOTONIC, &tiempo2);
   		diff_time(&tiempo2, &tiempo1, &latencia);
   		tiempoComputo = latencia.tv_sec*1000000000 + latencia.tv_nsec;  
		complejidad_temporal<< tiempoComputo<<" ";   
        operacion_23.push_back(complejidad_B);
        operacion_23.push_back(global_suma);
        operacion_23.push_back(global_resta);
        operacion_23.push_back(global_producto);
        Complejidad_23.insert(pair<int,vector<int> >(key,operacion_23));
        operacion_23.clear(); 
        for (int i =0;i<iterador->second.size(); i++){ // Se recorre hasta el final de los datos de las key.
            t.AgregarFreq(iterador->second[i],Mapa_Frecuencias[key][i]);
        }
        t.Recorrer_Arbol_23();// Se llama la funcion para guardar en un archivo el arbol 2-3 de datos.
        Archivo_23_datos<<endl;
        //Archivo_Complejidades23<<complejidad_B<<" "<<global_suma<<" "<<global_resta<<" "<<global_producto<<endl;
        Mapa_23_Datos.insert(pair<int, B_23  >(iterador->first,t.root)); // Se van ingresando en un map los arboles 2-3 creados.
        t.VolverNULLpuntero(); // Hacemos nulo el puntero del arbol 2-3.
        complejidad_B=0;
        global_suma=0;
        global_resta=0;
        global_producto=0;
        Archivo_23_datos<<endl;
        Archivo_23_datos<<endl;
        Archivo_23_freq<<endl;
        Archivo_23_freq<<endl;
    }
    cout<<"Arbol 2-3 Creado."<<endl;  
    complejidad_temporal<<endl;
	Archivo_23_datos.close();// Se cierra el archivo donde se escriben los arboles 2-3 de datos.
    Archivo_23_freq.close();// Se cierra el archivo donde guardamos los arboles 2-3 de frecuencias
//------------------------------------------------------------------------------------------------------------
	complejidad_B=0;
	complejidad_temporal<<"Arbol B"<<endl;
    Arbol_B BT(6);
	Archivo_b_datos.open("Arbol_B_Datos.dat"); // Se abre un archivo para guardar los arboles de datos B.
	Archivo_b_freq.open("Arbol_B_Frecuencias.dat"); // Se abre un archivo para guardar los arboles B de frecuencias de los datos.
    for(iterador =Mapa_Datos.begin(); iterador != Mapa_Datos.end(); ++iterador){// Se recorre el map de datos por key.
    	key=iterador->first;
    	Archivo_b_datos<<"Arbol B datos ()KEY: "<<key<<endl;
    	Archivo_b_freq<<"Arbol B frecuencias ()KEY: "<<key<<endl;
    	clock_gettime(CLOCK_MONOTONIC, &tiempo1);
        for (int i =0;i<iterador->second.size(); i++){// Hasta llegar al final de la linea de datos de cada key.
            BT.insert(iterador->second[i]); // Se insertan los datos a un arbol B.
        }
        clock_gettime(CLOCK_MONOTONIC, &tiempo2);
   		diff_time(&tiempo2, &tiempo1, &latencia);
   		tiempoComputo = latencia.tv_sec*1000000000 + latencia.tv_nsec;  
		complejidad_temporal<<  tiempoComputo<<" "; 
        operacion_b.push_back(complejidad_B);
        operacion_b.push_back(global_suma);
        operacion_b.push_back(global_resta);
        operacion_b.push_back(global_producto);
        Complejidad_b.insert(pair<int,vector<int> >(key,operacion_b));
        operacion_b.clear();		     
        for (int i =0;i<iterador->second.size(); i++){// Hasta llegar al final de la linea de datos de cada key.
            BT.AgregarFreq(iterador->second[i],Mapa_Frecuencias[key][i]);
        }        
        BT.Recorrer_Arbol_B();// Se llama a la funcion que escribe en un archivo los arboles B de datos creados.
        Archivo_b_datos<<endl;
        Archivo_b_datos<<endl;
        Archivo_b_freq<<endl;
        Archivo_b_freq<<endl;
        Mapa_BTREE_Datos.insert(pair<int, B_23  >(iterador->first,BT.root)); // Se van ingresando los arboles B creados a un map con sus key correspondientes.
        BT.VolverNULLpuntero(); // Volvemos nulo el puntero del arbol b.
        complejidad_B=0;
        global_suma=0;
        global_resta=0;
        global_producto=0;
    }
    cout<<"Arbol B Creado."<<endl; 
	cout<<endl; 
	Archivo_b_datos.close();// Se cierra el archivo donde se escribieron los arboles B de datos.
    Archivo_b_freq.close(); // Cerramos el archivo de escritura del arbol B de frecuencias.
    complejidad_temporal.close();
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
	int totales =0;
    map <int, ABB>:: iterator ti; // Se crea un iterador tipo ABB
    map <int, arbolAVL>:: iterator ti2; // Se crea un iterador tipo AVL
    map <int , B_23>:: iterator ti3;// Se crea un iterador tipo B y 2-3
    int r;
    long int Promedio=0;
    ofstream tiempos_ABB;
    tiempos_ABB.open("TiemposABB.dat"); // Se abre un archivo para guardar los tiempos de busqueda de los arboles ABB.
	ofstream tiempos_AVL;
    tiempos_AVL.open("TiemposAVL.dat");// Se abre un archivo para guardar los tiempos de busqueda de los arboles AVL.
    ofstream tiempos_23;
    tiempos_23.open("Tiempos23.dat");// Se abre un archivo para guardar los tiempos de busqueda de los arboles 2-3.
    ofstream tiempos_B;
    tiempos_B.open("TiemposB.dat");// Se abre un archivo para guardar los tiempos de busqueda de los arboles B.
    //Comienzo de las busquedas
//-----------------------------------------------------------------------------------------------------------------    
    //cout<<"Tiempo de arbol abb Datos: "<<endl;
    tiempos_ABB<<"Tiempo de arbol abb Datos: "<<endl;
    tiempos_ABB<<endl;
    long int PROM=0;
	Promedio =0;
    iterador = Busqueda.begin();
    for(ti = Mapa_ABB_Datos.begin(); ti != Mapa_ABB_Datos.end(); ++ti){ // Se recorre el mapa de aroboles abb de datos con un iterador.
        ABB_Aux=ti->second;// Se copia uno de los arboles en Prueba.
		key=ti->first;
		for (int i =0;i<iterador->second.size(); i++){
        	clock_gettime(CLOCK_MONOTONIC, &tiempo1); // Se toma el tiempo de inicio.
        	r=ABB_Tree_Busqueda(ABB_Aux, Busqueda[key][i]);// Se llama la funcion buscar con el arbol y el dato correspondiente ABB.
        	clock_gettime(CLOCK_MONOTONIC, &tiempo2);// Se toma el tiempo termino.
        	diff_time(&tiempo2, &tiempo1, &latencia); // Se saca la diferencia de ambos tiempos
        	tiempoComputo = latencia.tv_sec*1000000000 + latencia.tv_nsec; // Se consigue el tiempo en nanosegundos.
        	Promedio = Promedio + tiempoComputo;
        	PROM++;
        }
        Promedio = Promedio/PROM;
		++iterador;
        tiempos_ABB<<Promedio<<" "; // Se guarda dentro de un archivo los tiempos ABB.
		Promedio =0;
		PROM=0;
        totales++;
        ABB_Aux=NULL;
    }
    cout<<"Busquedas ABB Realizadas."<<endl;
    tiempos_ABB<<endl;
	iterador = Busqueda.begin();
//-----------------------------------------------------------------------------------------------------------------
	tiempos_ABB<<endl;
	tiempos_ABB<<"Final";
	tiempos_ABB.close();
	Promedio =0;
	PROM=0;
    tiempos_AVL<<"Tiempo de arbol AVL Datos: "<<endl;
    tiempos_AVL<<endl;
    for(ti2 = Mapa_AVL_Datos.begin(); ti2 != Mapa_AVL_Datos.end(); ++ti2){// Se recorre el map de arboles de datos AVL.
        AVL_Aux=ti2->second; // Se obtiene un arbol avl de datos
        key=ti2->first;
		for (int i =0;i<iterador->second.size(); i++){
        	clock_gettime(CLOCK_MONOTONIC, &tiempo1);// Se comienzan a tomar el tiempo de incio.
        	r=AVL_Tree_Busqueda(AVL_Aux,Busqueda[key][i]); // Mandamos a la funcion de busqueda el arbol y  el dato.
        	clock_gettime(CLOCK_MONOTONIC, &tiempo2);// Se saca el tiempo e termino.
        	diff_time(&tiempo2, &tiempo1, &latencia); // Se saca la diferencia de ambos tiempos para obtener el especifico.
        	tiempoComputo = latencia.tv_sec*1000000000 + latencia.tv_nsec; // Se obtiene el tiempo en nanosegundos.
        	Promedio = Promedio + tiempoComputo;
        	PROM++;
		}
        Promedio = Promedio/PROM;
		++iterador;
        tiempos_AVL<< Promedio<<" "; // Se guardan en un archivo los tiempos de buesquedas de datos avl.
 		Promedio =0;
		PROM=0;       
        totales++;
        AVL_Aux=NULL;
    }
    tiempos_AVL<<endl;
	tiempos_AVL<<endl;
	tiempos_AVL<<"Final";
	tiempos_AVL.close();
	cout<<"Busquedas AVL Realizadas."<<endl;
//-----------------------------------------------------------------------------------------------------------------
	iterador = Busqueda.begin();
 	Promedio =0;
	PROM=0; 	
    //cout<<"Tiempo de arbol 2-3 Datos: "<<endl;
    tiempos_23<<"Tiempo de arbol 2-3 Datos: "<<endl;
    tiempos_23<<endl;
    for(ti3 = Mapa_23_Datos.begin(); ti3 != Mapa_23_Datos.end(); ++ti3){ // Se recorre el mapa de arboles 2-3 de datos
        B_23_Aux=ti3->second; // Obtenemos un arbol 2-3.
		key=ti3->first;
		t.root=B_23_Aux;
		for (int i =0;i<iterador->second.size(); i++){
        	clock_gettime(CLOCK_MONOTONIC, &tiempo1); // Sacamos el tiempo de inicio.
        	t.search(Busqueda[key][i]);// Mandamos el numero a la busqueda de arbol 2-3 de datos.
        	clock_gettime(CLOCK_MONOTONIC, &tiempo2);// Tiempo de termino
        	diff_time(&tiempo2, &tiempo1, &latencia); // Diferencia entre tiempo de inicio y final.
        	tiempoComputo = latencia.tv_sec*1000000000 + latencia.tv_nsec; // obtencion de tiempo en nanosegundos
        	Promedio = Promedio+ tiempoComputo;
        	PROM++;
    	}
    	Promedio=Promedio/PROM;
        tiempos_23<< Promedio <<" "; // Se escribe en un archivo
        t.VolverNULLpuntero();
        B_23_Aux=NULL;
        Promedio=0;
        PROM=0;
        totales++;
        ++iterador;
    }
    tiempos_23<<endl;
	tiempos_23<<endl;
	tiempos_23<<"Final";
	tiempos_23.close(); // Se cierra el archivo 
	cout<<"Busquedas 2-3 Realizadas."<<endl;
//-----------------------------------------------------------------------------------------------------------------	
	iterador = Busqueda.begin();
    Promedio=0;
    PROM=0;	
    tiempos_B<<"Tiempo y complejidad de arbol b Datos: "<<endl;
    for(ti3 = Mapa_BTREE_Datos.begin(); ti3 != Mapa_BTREE_Datos.end(); ++ti3){// Se recorre el mapa de arboles b de datos.
        B_23_Aux=ti3->second; // Se obtiene un arbol b del map.
		key=ti3->first;
        BT.root=B_23_Aux;// Se iguala a una variable auxilair
        for (int i =0;i<iterador->second.size(); i++){
        	clock_gettime(CLOCK_MONOTONIC, &tiempo1);// se saca el tiempo de inicio.
        	BT.search(Busqueda[key][i]);
        	clock_gettime(CLOCK_MONOTONIC, &tiempo2);// Se saca el tiempo de termino.
        	diff_time(&tiempo2, &tiempo1, &latencia);// Se saca la diferencia de ambos tiempos.
        	tiempoComputo = latencia.tv_sec*1000000000 + latencia.tv_nsec; // Se obtiene el tiempo de busqueda en nanosegundos.
        	Promedio = Promedio + tiempoComputo;
        	PROM++;
    	}
    	Promedio = Promedio/PROM;
        tiempos_B<< Promedio << " "; // Se escriben en un rchivo los tiempos obtenidos.
        BT.VolverNULLpuntero();// Volvemos nulo en puntero del arbol B.
        Promedio =0;
        PROM=0;
        B_23_Aux=NULL;
        totales++;
        ++iterador;
    }
    cout<<"Busquedas B Realizadas."<<endl;
    tiempos_B<<endl;
	tiempos_B<<endl;
	tiempos_B<<"Final";
	tiempos_B.close();// Se cierra el archivo donde se guardan los tiempos de busqueda.
//-----------------------------------------------------------------------------------------------------------------
	cout<<totales<<endl;

}

int abrir_archivo()
{
    cout<<"Ejecutando el programa..."<<endl;
    cout<<endl;
	vector<int> Datos;
	vector<int> Vector_Principal; // Vector de datos sin datos repetidos
	vector<int> Vector_Frecuencia;// Vector con las frecuencias de los datos del vector de arriba.
	map <int , vector<int> > Mapa_Datos; // Mapa de los datos del archivo.
	map <int , vector<int> > Map_Datos;// Mapa de los datos del archivo.
	map <int , vector<int> > Mapa_Frecuencias;// Mapa de las frecuencias de los datos del archvo.
	map <int , vector<int> > Busqueda;
	ifstream Archivo; // funcion abrir archivo.
	Archivo.open("ubLista3Mb_Factor2_64.dat"); // Se abre el archivo con miles de datos.
	int key; // se declara la key de los maps
	int dato;
	int DATOS_TOTALES; // cantidad de datos en una fila.
	int contador=0;	
	int inservible; // dato inservible dentro de la linea.
	int Numero_Final_Busqueda = 0;
	int numero;
	while(!Archivo.eof()) // mientras la lectura del archivo no llegue al final.
	{	
		Archivo>>key; // Recibe la key que es el primer elemento de una fila.
		Archivo>>inservible; // recibe el segundo elemento de la fila que no sirve.
		Archivo>>inservible;// recibe el tercer elemento de la fila que tampoco sirve.
		Archivo>>DATOS_TOTALES; // recibe el cuarto elemento de la fila que es igual a la cantidad de datos de la fil.
		contador=0;
		while(DATOS_TOTALES!=contador) // Mietras contador sea distinto a la cantidad de datos de una fila.
		{
			Archivo>>dato; // Se lee un dato de la fila
			Datos.push_back(dato); // se guarda en un vector de datos.
			contador++; // el contador aumenta en 1.
		}
		if(DATOS_TOTALES!=0){
			Busqueda.insert(pair<int,vector<int> >(key,Datos));
			Map_Datos.insert(pair<int,vector<int> >(key,Datos));	
			for(int i=0; i<DATOS_TOTALES;i++){ // for encargado de crear un vector de datos no repetidos
				if(i==1)
				{
					Vector_Principal.push_back(Datos[0]); // agrega al vector pincipal el primer elemento de la fila.
					Numero_Final_Busqueda++;
				}
				else{
					numero = Datos[i];
					numero=esta(numero,Vector_Principal,Numero_Final_Busqueda);// Se envia a una fncion que retorna si esta el elemento repetido.
					if(numero==0){
						Vector_Principal.push_back(Datos[i]);// Asi se crea un vector con elementos no repetidos para poder sacar las frecuencias.
						Numero_Final_Busqueda++;					
					}
				}
			}
			for(int i=0; i<Numero_Final_Busqueda; i++) // For para sacar las frecuencias de los datos del vector.
			{
				numero=Vector_Principal[i];
				numero=frecuencia(numero,Datos, DATOS_TOTALES); // Funcion que retorna la frecuencia de un dato de nuestro vector sin repeticiones.
				Vector_Frecuencia.push_back(numero); // Se agrega a la misma posicion de otro vector el dato y su frecuencia.
			}
			Mapa_Datos.insert(pair<int,vector<int> >(key,Vector_Principal)); // se agrega al map de datos su key y su vector de datos.
			Mapa_Frecuencias.insert(pair<int,vector<int> >(key,Vector_Frecuencia));// Se crea un map con la key aterior pero con un vector de frecuencias.
			Numero_Final_Busqueda=0;
			Datos.clear();// Se limpan los vectores.
			Vector_Principal.clear();// Se limpan los vectores.
			Vector_Frecuencia.clear();// Se limpan los vectores.
		}
	}
	cout<<endl;
	cout<<endl;
	cout<<"Archivo cargado Completamente"<<endl;
	cout<<endl;	
	Crear_Arboles(Mapa_Datos,Mapa_Frecuencias, Busqueda, Map_Datos);// Se llama a la funcion para la creacion de arboles.
    cout<<"Arboles creados Exitosamente"<<endl;
	Archivo.close(); // Se cierra el archivo de los miles de datos.
	return 0;
}

int main(){ // Funion main principal.
	abrir_archivo();
	return 0;
}
