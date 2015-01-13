#ifndef _STRUCTURE_
#define _STRUCTURE_

#define MAX_HEADER_LENGTH (500)

class t_structure
{
public:
    int numofbases;
	int* numseq;
	char* nucs;
	int* basepr;
	char* ctlabel;
	char* fp;

	t_structure(char* fp); // Can take a seq or a ct file.
	t_structure(); // Defaults constructor.
	~t_structure();

	void openct(char* ct_fp);
	void openseq(char* seq_fp);

	bool verify_ct(char* ct_fp);
	bool verify_seq(char* seq_fp);
};

#endif // _STRUCTURE_

