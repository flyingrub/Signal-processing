#ifndef WAVE_HPP
#define WAVE_HPP

class Wave {


private :
  //chunk type
  char file_type[4];         // (4 octets) : Constante "RIFF" i.e identification du format
  int file_size;             // (4 octets) : file_size est le nombre d'octet restant � lire (i.e = taille du fichier moins 8 octets)
  char file_id[4];           // (4 octets) : Identifiant "WAVE"
  //chunk format
  char chunk_id[4];          // (4 octets) : Identifiant "fmt "
  int chunk_size;            // (4 octets) : Nombre d'octets utilis�s pour d�finir en d�tail le chunk
  short format;              // (2 octets) : Format de fichier (1: PCM,  ...)
  short channels_nb;         // (2 octets) : Nombre de canaux (1 pour mono ou 2 pour st�r�o)
  int sampling_freq;         // (4 octets) : Fr�quence d'�chantillonnage (en Hertz)
  int bytes_per_second;      // (4 octets) : Nombre d'octets par seconde de musique
  short bytes_per_sample;    // (2 octets) : Nombre d'octets par �chantillon
  short depth;               // (2 octets) : Nombre de bits par donn�e (8 ou 16)
  //chunk donn�es
  char data_id[4];          // (4 octets) : Constante "data"
  int data_size;            // (4 octets) : nombre d'octet restant (i.e = taille du fichier moins 44 octets)
  //data
  unsigned char* data8;     // Tableau de don�es lorsque l'on est sur des donn�es 8 bits
  short* data16;            // Tableau de don�es lorsque l'on est sur des donn�es 16 bits
  long int data_nb;         // Nombre de donn�es

  bool is_data8_allocated;  // Indique si le tableau � �t� allou�
  bool is_data16_allocated; // Indique si le tableau � �t� allou�

  static void checkTypesSize(); // Test la taille des types

  void InitDescriptor(int depth,          // Nombre de bits par donn�e (8 ou 16)
                      short channels_nb,  // Nombre de canaux (1 pour mono ou 2 pour st�r�o)
                      int sampling_freq); // Fr�quence d'�chantillonnage (en Hertz)(classique en musique : 44100 Hz)

public:

  //M�thode permettant de r�cup�rer les donn�es (cad le signal sur 8 bits)
  void getData8(unsigned char** data, // pointeur sur le tableau de don�es lorsque l'on est sur des donn�es 8 bits
                int* size);           // pointeur sur la taille du tableau
  //M�thode permettant de modifier les donn�es
  void modifData8(unsigned char* data); // Tableau de don�es lorsque l'on est sur des donn�es 8 bits

  void read(char* fileName);   //Lecture d'un fichier wave
  void write(char* fileName);

  int getSampleRate();

  Wave();
  ~Wave();
  Wave(short* data16,      // Tableau de donn�es
       long int data_nb,   // Nombre de donn�es
       short channels_nb,  // Nombre de canaux (1 pour mono ou 2 pour st�r�o)
       int sampling_freq); // Fr�quence d'�chantillonnage (en Hertz)(classique en musique : 44100 Hz)
  Wave(unsigned char* data8, // Tableau de don�es lorsque l'on est sur des donn�es 8 bits
       long int data_nb,     // Nombre de donn�es
       short channels_nb,    // Nombre de canaux (1 pour mono ou 2 pour st�r�o)
       int sampling_freq);   // Fr�quence d'�chantillonnage (en Hertz)(classique en musique : 44100 Hz)
};
#endif //IO_WAVE_HPP
