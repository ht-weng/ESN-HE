#include <iostream>

#include <helib/FHE.h>
#include <helib/EncryptedArray.h>

int main(int argc, char *argv[]) {
  /*  Example of BGV scheme  */
  
  // Plaintext prime modulus
  unsigned long p = 4999;
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = 32109;
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of bits of the modulus chain
  unsigned long bits = 300;
  // Number of columns of Key-Switching matrix (default = 2 or 3)
  unsigned long c = 2;
  
  std::cout << "Initialising context object..." << std::endl;
  // Intialise context
  FHEcontext context(m, p, r);
  // Modify the context, adding primes to the modulus chain
  std::cout  << "Building modulus chain..." << std::endl;
  buildModChain(context, bits, c);

  // Print the context
  context.zMStar.printout();
  std::cout << std::endl;
  
  // Print the security level
  std::cout << "Security: " << context.securityLevel() << std::endl;
  
  // Secret key management
  std::cout << "Creating secret key..." << std::endl;
  // Create a secret key associated with the context
  FHESecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();
  std::cout << "Generating key-switching matrices..." << std::endl;
  // Compute key-switching matrices that we need
  addSome1DMatrices(secret_key);
  
  // Public key management
  // Set the secret key (upcast: FHESecKey is a subclass of FHEPubKey)
  const FHEPubKey& public_key = secret_key;
  
  // Get the EncryptedArray of the context
  const EncryptedArray& ea = *(context.ea);
  
  // Get the number of slot (phi(m))
  long nslots = ea.size();
  std::cout << "Number of slots: " << nslots << std::endl;
  
  // Create a vector of long with nslots elements
  std::vector<long> ptxt(nslots);
  // Set it with numbers 0..nslots - 1
  for (int i = 0; i < nslots; ++i) {
    ptxt[i] = i;
  }

  // PlaintextArray p0(ea);
  // std::vector<double> p0(2);
  double p0 [] = { 16, 2, 77, 40, 12071 };
  // p0 = [0.0 0.1];

  // Print the plaintext
  std::cout << "Initial Ptxt: " << ptxt << std::endl;
  std::cout << "Initial Ptxt: " << p0 << std::endl;
  
  // Create a ciphertext
  Ctxt ctxt(public_key);
  Ctxt c0(public_key);
  // Encrypt the plaintext using the public_key
  ea.encrypt(ctxt, public_key, ptxt);
  ea.encrypt(c0, public_key, p0);

  // Square the ciphertext
  // ctxt *= ctxt;
  // Double it (using additions)
  // ctxt += ctxt;
  
  // Create a plaintext for decryption
  std::vector<long> decrypted(nslots);
  // PlaintextArray pp0(ea);
  double pp0 [] = { 2, 2, 77, 40, 12071 };
  // Decrypt the modified ciphertext
  ea.decrypt(ctxt, secret_key, decrypted);
  ea.decrypt(c0, secret_key, pp0);
  
  // Print the decrypted plaintext
  std::cout << "Decrypted Ptxt: " << decrypted << std::endl;
  std::cout << "Decrypted Ptxt: " << pp0 << std::endl;
  
  return 0;
}
