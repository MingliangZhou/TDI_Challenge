#include <iostream>
#include <cmath>

using namespace std;

class Card
{
  public:
    Card(int ncard = 1, int nsuit = 1);
    virtual ~Card() { reset(); }

    void set_ptable(int ncard, int nsuit);

    long double get_cdf(int points);
    long double get_mean();
    long double get_var();

  protected:
    void reset();
    void get_ncut(int index, int ncut[]);
    void fill_factorial();
    void fill_ncomb();
    void fill_nsuitcut();
    void fill_ptable();

    int N;
    int M;
    int K;
    int KM;

    long double *factorial;
    long double **ncomb;
    long double *nsuitcut;
    long double *ptable;
};

Card::Card(int ncard, int nsuit):
  factorial(NULL),
  ncomb(NULL),
  nsuitcut(NULL),
  ptable(NULL)
{
  set_ptable(ncard, nsuit);
}

void Card::set_ptable(int ncard, int nsuit)
{
  if(ncard <= 0 || nsuit <= 0 || ncard % nsuit != 0)
  {
    cerr << "Wrong input for (N,M) = (" << ncard << "," << nsuit <<")" << endl;
    return;
  }

  reset();

  N = ncard;
  M = nsuit;
  K = N / M;
  KM = pow(K,M);

  factorial = new long double[N+M+1];
  ncomb = new long double*[N];
  for(int i=0; i<N; i++)
    ncomb[i] = new long double[N];
  nsuitcut = new long double[KM];
  ptable = new long double[N];

  fill_factorial();
  fill_ncomb();
  fill_nsuitcut();
  fill_ptable();

  return;
}

long double Card::get_cdf(int points)
{
  if(points >= N-M) return 1.;

  long double psum = 0.;
  for(int i=0; i<=points; i++)
    psum += ptable[i];

  return psum;
}

long double Card::get_mean()
{
  long double mean = 0.;
  for(int i=0; i<=N-M; i++)
    mean += i * ptable[i];

  return mean;
}

long double Card::get_var()
{
  long double mean = get_mean();
  long double var = 0.;
  for(int i=0; i<=N-M; i++)
    var += i*i * ptable[i];
  var -= mean*mean;

  return var;
}

void Card::reset()
{
  if(factorial)
    delete[] factorial;
  for(int i=0; i<N; i++)
    if(ncomb && ncomb[i])
      delete[] ncomb[i];
  if(nsuitcut)
    delete[] nsuitcut;
  if(ptable)
    delete[] ptable;

  N = 1;
  M = 1;
  K = 1;
  KM = 1;

  return;
}

void Card::get_ncut(int index, int ncut[])
{
  for(int i=0; i<M; i++)
    ncut[i] = 0.;

  int i = 0;
  while(index > 0)
  {
    ncut[i] = index % K;
    index /= K;
    i++;
  }

  return;
}

void Card::fill_factorial()
{
  for(int i=0; i<=N+M; i++)
  {
    factorial[i] = 1.;
    for(int j=2; j<=i; j++)
      factorial[i] *= j;
  }

  return;
}

void Card::fill_ncomb()
{
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
    {
      if(j <= i)
        ncomb[i][j] = factorial[i] / factorial[j] / factorial[i-j];
      else
        ncomb[i][j] = 0.;
    }

  return;
}

void Card::fill_nsuitcut()
{
  int *ncut = new int[M];
  int *lcut = new int[M];

  for(int i=0; i<KM; i++)
  {
    get_ncut(i, ncut);

    bool sym = false;
    for(int j=0; j<M; j++)
    {
      if(sym) break;
      for(int k=0; k<j; k++)
        if(ncut[k] < ncut[j])
        {
          int li = i + (ncut[j]-ncut[k])*pow(K,k) + (ncut[k]-ncut[j])*pow(K,j);
          nsuitcut[i] = nsuitcut[li];
          sym = true;
          break;
        }
    }
    if(sym) continue;

    int ncut_sum = 0;
    for(int j=0; j<M; j++)
      ncut_sum += ncut[j];

    long double ncut_prod = 1.;
    for(int j=0; j<M; j++)
      ncut_prod *= factorial[ncut[j]+1];

    nsuitcut[i] = factorial[ncut_sum+M] / ncut_prod;
    for(int j=0; j<i; j++)
    {
      get_ncut(j, lcut);

      long double lcut_prod = 1.;
      for(int k=0; k<M; k++)
        lcut_prod *= ncomb[ncut[k]][lcut[k]];

      nsuitcut[i] -= lcut_prod * nsuitcut[j];
    }
  }

  delete[] ncut;
  delete[] lcut;
  return;
}

void Card::fill_ptable()
{
  int *ncut = new int[M];

  for(int i=0; i<KM; i++)
  {
    get_ncut(i, ncut);

    int ncut_sum = 0;
    for(int j=0; j<M; j++)
      ncut_sum += ncut[j];

    long double ncut_prod = 1.;
    for(int j=0; j<M; j++)
      ncut_prod *= ncomb[K-1][ncut[j]];

    int points = N - M - ncut_sum;
    ptable[points] += nsuitcut[i] * ncut_prod;
  }

  const long double ntotal = factorial[N] / pow(factorial[K],M);
  long double psum = 0.;
  for(int i=0; i<N; i++)
  {
    ptable[i] /= ntotal;
    psum += ptable[i];
    cout << i << " points, prob = " << ptable[i] << endl;
  }
  cout << "Sum of probability for (N,M) = (" << N << "," << M <<") is " << psum << endl;

  delete[] ncut;
  return;
}

int main()
{
  cout.precision(15);

  int N = 26, M = 2;
  Card *card = new Card(N,M);
  cout << "Mean for (N,M) = (" << N << "," << M <<") is " << card->get_mean() << endl;
  cout << "Standard deviation for (N,M) = (" << N << "," << M <<") is " << sqrt(card->get_var()) << endl;
  cout << "Conditional probability (points>12|points>6) for (N,M) = (" << N << "," << M <<") is " << (1.-card->get_cdf(12))/(1.-card->get_cdf(6)) << endl;

  N = 52, M = 4;
  card->set_ptable(N,M);
  cout << "Mean for (N,M) = (" << N << "," << M <<") is " << card->get_mean() << endl;
  cout << "Standard deviation for (N,M) = (" << N << "," << M <<") is " << sqrt(card->get_var()) << endl;
  cout << "Conditional probability (points>12|points>6) for (N,M) = (" << N << "," << M <<") is " << (1.-card->get_cdf(12))/(1.-card->get_cdf(6)) << endl;

  delete card;
  return 0;
}
