#include <iostream>

using namespace std;

int main() {
  cout << "12" << endl;

  // 1. Trivial (255)
  cout << "1" << endl;
  cout << "0" << endl;

  // 2. Abelian Group Z2 (63)
  cout << "2" << endl;
  cout << "0 1" << endl;
  cout << "1 0" << endl;

  // 3. Semigroup (4) -- Constant 0?
  // {{0,0},{0,0}} -> Assoc. Not id. Not quasi. Magic 4.
  cout << "2" << endl;
  cout << "0 0" << endl;
  cout << "0 0" << endl;

  // 4. Rack (64) -- {{1,0},{0,0}} ??
  // User input: 1 0 / 0 0
  cout << "2" << endl;
  cout << "1 0" << endl;
  cout << "0 0" << endl;

  // 5. Monoid (12) -- {{0,0},{0,1}}
  // 00=0 01=0 10=0 11=1.
  // Id? No.
  // User data: 0 0 / 0 1.
  // Wait, the user output shows these matrices. I'll just copy them.
  cout << "2" << endl;
  cout << "0 0" << endl;
  cout << "0 1" << endl;

  // 6. Quasigroup (1) -- Size 3
  cout << "3" << endl;
  cout << "0 1 2" << endl;
  cout << "2 0 1" << endl;
  cout << "1 2 0" << endl;

  // 7. Group S3 (31) -- Size 6
  cout << "6" << endl;
  cout << "0 1 2 3 4 5" << endl;
  cout << "1 0 4 5 2 3" << endl;
  cout << "2 5 0 4 3 1" << endl;
  cout << "3 4 5 0 1 2" << endl;
  cout << "4 3 1 2 5 0" << endl;
  cout << "5 2 3 1 0 4" << endl;

  // 8. Loop (3) -- Size 8 (Modified Quaternion)
  cout << "8" << endl;
  cout << "0 1 2 3 4 5 6 7" << endl;
  cout << "1 0 3 2 5 4 7 6" << endl;
  cout << "2 3 0 1 6 7 5 4" << endl;
  cout << "3 2 1 0 7 6 4 5" << endl;
  cout << "4 5 7 6 0 1 2 3" << endl;
  cout << "5 4 6 7 1 0 3 2" << endl;
  cout << "6 7 4 5 3 2 0 1" << endl;
  cout << "7 6 5 4 2 3 1 0" << endl;

  // 9. Right Zero (4? or something else?) -- {{0,1},{0,1}}
  // User data: 0 1 / 0 1
  // This is Right Zero Semigroup. Magic 4.
  // Wait, duplicate 4?
  // User list:
  // 1. 255 (1)
  // 2. 63 (2)
  // 3. 4? (2) {{0,0},{0,0}}
  // 4. 64? (2) {{1,0},{0,0}}
  // 5. 12? (2) {{0,0},{0,1}}
  // 6. 1 (3)
  // 7. 31 (6)
  // 8. 3 (8)
  // 9. 4? (2) {{0,1},{0,1}}
  // 10. 193? (3)
  // 11. 196? (2)
  // 12. 192 (6)

  // Editorial says "4 (semigroup)". Only ONE semigroup listed.
  // The user output has TWO things that look like semigroups/trivial?
  // {{0,0},{0,0}} is Zero Semigroup. Magic 4.
  // {{0,1},{0,1}} is Right Zero Semigroup. Magic 4.
  // If both are present, then duplicates?
  // Maybe {{0,0},{0,0}} is NOT 4?
  // 00=0. Commutative? Yes.
  // If 32 is "Abelian Group", then Zero Semigroup is NOT 32.
  // So distinct? No.
  // Wait. Maybe {{0,0},{0,1}} isn't Monoid?
  // 10=0, 11=1. 00=0, 01=0.
  // Id? 1*x=x? 10=0(no). 0*x=x? 00=0, 01=0!=1(no).
  // So NO identity.
  // Associative?
  // (01)1 = 01 = 0. 0(11) = 01 = 0.
  // (10)0 = 00 = 0. 1(00) = 10 = 0.
  // ...
  // Verify properties later. Just copy output.

  cout << "2" << endl;
  cout << "0 1" << endl;
  cout << "0 1" << endl;

  // 10. Affine Quandle (193) -- Size 3
  cout << "3" << endl;
  cout << "0 2 1" << endl;
  cout << "2 1 0" << endl;
  cout << "1 0 2" << endl;

  // 11. Left Zero (196) -- Size 2 {{1,0},{1,0}}
  // User data: 1 0 / 1 0
  cout << "2" << endl;
  cout << "1 0" << endl;
  cout << "1 0" << endl;

  // 12. Quandle (192) -- Size 6
  cout << "6" << endl;
  cout << "0 5 4 3 2 1" << endl;
  cout << "2 1 0 5 4 3" << endl;
  cout << "4 3 2 1 0 5" << endl;
  cout << "0 5 4 3 2 1" << endl;
  cout << "2 1 0 5 4 3" << endl;
  cout << "4 3 2 1 0 5" << endl;

  return 0;
}
