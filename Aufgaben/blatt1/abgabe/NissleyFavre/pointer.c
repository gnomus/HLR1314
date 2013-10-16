#include <stdio.h>

int TODO;

void basic_pointer(int x)
{
	int *adresse_von_x;

	adresse_von_x = &x;

	printf("Der Wert von x ist: %d\n", x /* TODO */);
	printf("Die Adresse von x ist %p\n", &x /* TODO */);
	printf("Adresse von x mittels adresse_von_x %p\n", adresse_von_x /* TODO */);
	printf("Wert von x mittels adresse_von_x: %d\n", *adresse_von_x /* TODO */);

}

void basic_pointer2(int x)
{
	int *adresse_von_x;
	adresse_von_x = &x /* TODO */;

	//Eine andere Variable y erhaelt den Wert von x
	int y = x /* TODO */;

	printf("Der Wert von y ist %d\n", y /* TODO */);

	//Zuweisung über Adresse
	x = 10;
	y = *adresse_von_x;
	printf("Der Wert von y ist %d\n", y /* TODO */);
}

void basic_pointer_changeValue(int x)
{
	int *adresse_von_x;
	adresse_von_x = &x /* TODO */;
	//Ändern Sie den Wert von x zu 10
	x = 10;
	printf("x = %d\n", x /* TODO */);

	//Ändern Sie den Wert von x über seine ADRESSE
	*adresse_von_x = 20;
	printf("x = %d\n", x);

}


void call_by_reference(int *x)
{
	//Ändern Sie den Wert, der an der Adresse steht, die im Wert x gespeichert ist.
	*x = 200;
}

int main()
{
	int x = 5;

	basic_pointer(x);
	basic_pointer2(x);
	basic_pointer_changeValue(x);
	printf("Wert von x vor der Funktion call_by_reference: %d\n",x);
	call_by_reference(&x);
	printf("Wert von x nach der Funktion call_by_reference: %d\n",x);

	return 0;

}
