#include <iostream>
/*
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/EcuyerCombined.h"
#include "../RandomGenerator/FiniteSet.h"
#include "../RandomGenerator/Exponential.h"
#include "../RandomGenerator/Normal.h"
#include "../RandomGenerator/Poisson.h"
#include "../SDE/BSEuler1D.h"
#include "../PDE/PDEGrid2D.h"
*/
void TestPDE();
void TestRandom();
void TestSDE();

int main()
{
    //TestRandom();
    //TestSDE();
    //TestPDE();    
	std::cout << "hello world" << std::endl;
}

void TestPDE()
{
    /*
        Evaluate a Call option in the Black-Scholes framework:
        Expected Price : 6.80495771
    */
    double Spot = 100;
    double Strike = 100;
    double Maturity = 1;
    double Rate = 0.05;
    double Volatility = 0.1;

    double SMin = 0;
    double SMax = 1000;

    R2R1Function* VarianceFunction = new BSVariance(Volatility);
    R2R1Function* TrendFunction = new BSTrend(Rate);
    R2R1Function* ActualizationFunction = new BSActualization(Rate);
    R2R1Function* SourceTermFunction = new NullFunction();

    R1R1Function* TopBoundaryFunction = new CallTopBoundary(SMax, Strike);
    R1R1Function* BottomBoundaryFunction = new CallBottomBoundary(SMin, Strike);
    R1R1Function* RightBoundaryFunction = new CallTerminalCondition(Strike);

    PDEGrid2DExplicit BlackScholesGrid(Maturity, SMin, SMax, 10000, 1.,
        VarianceFunction, TrendFunction, ActualizationFunction, SourceTermFunction,
        TopBoundaryFunction, BottomBoundaryFunction, RightBoundaryFunction);

    BlackScholesGrid.FillNodes();
    double priceAtZero = BlackScholesGrid.GetTimeZeroNodeValue(Spot);

    std::cout << "Price = " << priceAtZero << std::endl;
}

void TestSDE()
{
    myLong Multiplier, Increment, Modulus, Seed;

    Multiplier = 40014;
    Increment = 0;
    Modulus = 2147483563;
    Seed = 1;

    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);

    BSEuler1D Scheme = BSEuler1D(NormBox, 100., 0.5, 0.1);
    Scheme.Simulate(0., 1., 100);
    SinglePath* path = Scheme.GetPath(0);
    vector<double>& pathValues = path->GetValues();
    vector<double>& pathTimes = path->GetTimes();

    for (size_t i = 0; i < pathTimes.size(); ++i)
    {
        std::cout << pathTimes[i] << " " << pathValues[i] << std::endl;
    }

    /*
    PoissonFirstAlgo* PoissonFirst = new PoissonFirstAlgo(0.5, Unif);
    isAverageOk = PoissonFirst->TestMean(100000, 0.01);
    isVarianceOk = PoissonFirst->TestVariance(100000, 0.1);
    */
    // isAverageOk = X->TestMean(100000, 0.01);
    // isVarianceOk = X->TestVariance(100000, 0.01);

    /*
    bool isOk = X.TestMean(100000, 0.000001);
    std::cout << ((isOk) ? "OK" : "Not good") << std::endl;
    */
}

void TestRandom()
{
    myLong Multiplier, Increment, Modulus, Seed;

    Multiplier = 40014;
    Increment = 0;
    Modulus = 2147483563;
    Seed = 1;

    //X = new LinearCongruential(Multiplier, Increment, Modulus, Seed);
    //X = new EcuyerCombined();
    UniformGenerator* Unif = new EcuyerCombined();
    std::vector<double> probs(3);
    probs[0] = 0.4;
    probs[1] = 0.5;
    probs[2] = 0.1;
    RandomGenerator* Finit = new FiniteSet(probs, Unif);
    /*
        for (size_t i = 0; i < 10; ++i)
        {
            std::cout << Finit->Generate() << std::endl;
        }

        ExponentialInverseDistribution* Exp = new ExponentialInverseDistribution(0.5, Unif);
        */
    bool isAverageOk, isVarianceOk;
    //isAverageOk = Exp->TestMean(100000, 0.01);
    //isVarianceOk = Exp->TestVariance(100000, 0.1);

    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif);
    isAverageOk = NormBox->TestMean(100000, 0.01);
    isVarianceOk = NormBox->TestVariance(100000, 0.1);
}