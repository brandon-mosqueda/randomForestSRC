#include  <stdlib.h>
#include  <math.h>
#include  "splitCustom.h"

#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
#include <stdio.h>

/*
  Before coding a custom split rule, you must register it with the
  family specific vector of custom split rules.  This allows the user
  to code from one (1) to sixteen (16) different split rules for each
  family.  In the following function, each sample split rule is
  registered in slot one (1) of the family specific vector of custom
  split rules.

  The specific rule for slot X in the vector, where 1 <= X <= 16, is
  activated by the R-code user parameter:  splitrule = "customX".

  Note that splitrule = "custom" is the equivalent of splitrule = "custom1".
 */

void registerCustomFunctions()
{
	// Register the custom classification split rule in the first slot.
	registerThis (&getCustomSplitStatisticMultivariateClassification,
				  CLAS_FAM, 1);

	// Register the custom regression split rule in the first slot.
	//registerThis (&cost, REGR_FAM, 1);
	registerThis (&getCustomSplitZipPoisson, REGR_FAM, 1);
	registerThis (&getCustomSplitStatisticMultivariateRegression, REGR_FAM, 2);

	// Register the custom survival split rule in the first slot.
	registerThis (&getCustomSplitStatisticSurvival, SURV_FAM, 1);

	// Register the custom competing risk split rule in the first slot.
	registerThis (&getCustomSplitNonHomogeneous, CRSK_FAM, 1);
	registerThis (&getCustomSplitStatisticCompetingRisk, CRSK_FAM, 2);


	// If you have more than one function of each type, you would
	// uncomment one or more of the following statements and complete
	// the coding for the named functions.  Note that you can call these
	// functions anything you want.  In a multivariate scenario you must
	// be sure that the index of both the classification and regression
	// split rule are identical.  For example, in a multivariate
	// scenario, splitrule = "custom2", will call the classification
	// split rule in the second slot, along with the regression split
	// rule in the second slot if the respones contain factors and reals.

	//  registerThis (&getCustomSplitStatisticMultivariateClassificationTwo,
	//  			  CLAS_FAM, 2);
	//  registerThis (&getCustomSplitStatisticMultivariateRegressionTwo,
	//  			  REGR_FAM, 2);
	//  registerThis (&getCustomSplitStatisticSurvivalTwo, SURV_FAM, 2);
	//  registerThis (&getCustomSplitStatisticCompetingRiskTwo, CRSK_FAM, 2);
}

/*
  Generic Custom Split Rule Harness

  FUNCTION INPUTS:
  n - number of replicates in the parent node to be split

  membership     - vector of length n representing daughter node
  membership (LEFT or RIGHT) as a result of the
  test split.

  time     - vector of length n of time variable for survival data sets.
  This will be NULL for non-survival data sets.
  event   - vector of length n of event/censoring variable for
  survival data sets.  This will be NULL for non-survival
  data sets.

  response - vector of length n of response variable (y-outcome) in the
  node. In the multivariate case this function will be called
  once for each response.  The multivariate split statistic
  is the sum of the individual response statistic.  It is
  important to normalize the statistic for each response to
  avoid one response overwhelming another's value.  This will
  be NULL for survival data sets.

  mean     - convenience value representing the mean of the response vector
  variance - convenience value representing the variance of the response vector

  maxLevel - convenience value representing the maximum level for this factor
  respone in the data set as a whole.  This will be NULL for non-factor
  respones.

  feature - matrix of user specified features that can be sent into
  the split rule and acted on as desired.  The matrix is of dimension
  [featureCount] x [n].  Features are neither y-variables or
  x-variables. However, for expediency, they are specified by the user
  as y-variables, but are tagged as having zero weight via the
  y-variable weight vector.  Thus, they are never used in the the
  pre-defined split rules, and have no predicted value. The pointer will be NULL
  when features are absent.

  featureCount - count of features in the above matrix, specifically the
  number of rows.  The count will be zero when features are absent.


  FUNCTION OUTPUT:  returns a double value representing the custom split
  statistic

  Depending on the split statistic, the user may need to allocate and
  de-allocate arrays of various dimensions.  A typical alloc/de-alloc is
  defined in this file:

  unsigned int *alloc_uivector(...)
  void        dealloc_uivector(...)

  Always remember to deallocate what has been allocated.

  Always remember to declare the function in the corresponding ".h" file.
*/

// Newton Raphson Function to resolve a non-linear equation
double new_raph(int Nplus, int ntotal)
{
	// 0.0001 is the error level we wish
    double tol = 0.000001;
    int iter = 0;
    double x_old, x_new;

    if (ntotal == 0)
    {
    	x_old = 0.;
    }
    else
    {
		//x_old = 1 / (double)(Nplus);
    	x_old = -log((double)(Nplus) / (double)(ntotal)) /
    			(double)(ntotal) + 0.2;
    	do
    	{
    		x_new = x_old;
    		x_old = x_new - (x_new * (double)(Nplus) +
    						 exp(-x_new * (double)(ntotal)) - 1) /
    						((double)(Nplus) - (double)(ntotal) *
    						 exp(-x_new*(double)(ntotal)));
    		iter += 1;
    	} while(fabs(x_old - x_new) > tol && iter < 100);

    	return x_old;
    }
}

double getCustomSplitZipPoisson(
	unsigned int n,
	char *membership,
	double *time,
	double *event,

    unsigned int eventTypeSize,
    unsigned int eventTimeSize,
    double *eventTime,

    double *response,
    double mean,
    double variance,
    unsigned int maxLevel,

	double **feature,
	unsigned int featureCount)
{
	// EXAMPLE:  Multivariate Regression

	// Local variables needed for this example:
	double sumLeft, sumRght;
	double sumLeftSqr, sumRghtSqr;
	double delta;

	// Left and right normalization sizes.
	unsigned int leftSize, rghtSize;

	unsigned int i;

	// Initialization of local variables:
	sumLeft = sumRght = 0.0;
	leftSize = rghtSize = 0;

	delta = 0.0;

	int sum_eve_left;
	sum_eve_left = 0;

	int sum_eve_right;
	sum_eve_right = 0;

	int N0_left, Nplus_left, N0_right, Nplus_right;
	N0_left = Nplus_left = N0_right = Nplus_right = 0;

	double alpha_t_left;
	alpha_t_left = 0;

	double alpha_t_right;
	alpha_t_right = 0;

	double const_left;
	double const_right;

	// In general, calculating a split statistic will require iterating
	// over all members in the parent node, and ascertaining daughter
	// membership, and performing a well defined calculation based on
	// membership.  In this example, the sum of the difference from the
	// mean for the y-outcome in each daughter node is calculated.

	for (i = 1; i <= n; i++)
	{
		// Memconst_leftconst_leftbership will be either LEFT or RIGHT.
		if (membership[i] == LEFT)
		{
			// Add the left member to the sum.
			sum_eve_left = response[i];

			if(sum_eve_left == 0)
				N0_left += 1;
			else
				Nplus_left += 1;

			sumLeft += response[i];
			leftSize++;
		}
		else
		{
			sum_eve_right = response[i];
			if (sum_eve_right == 0)
				N0_right += 1;
			else
				Nplus_right += 1;

			sumRght += response[i];
			rghtSize++;
		}
	}

	const_left = new_raph(Nplus_left, sumLeft);
	alpha_t_left = (double)(sumLeft) * const_left;

	const_right = new_raph(Nplus_right, sumRght);
	alpha_t_right = (double)(sumRght) * const_right;

	// Finally, we calculate the composite mean square error for each daughter.

	if(alpha_t_left > 0.000001)
		sumLeftSqr += (double)(sumLeft) * log(alpha_t_left);
	sumLeftSqr += -(double)(Nplus_left) * alpha_t_left;

	if(!sumLeft == 0)
		sumLeftSqr += -(double)(Nplus_left) * log(1. - exp(-(double)(sumLeft) *
														   const_left));

	if(alpha_t_right > 0.000001)
		sumRghtSqr += (double)(sumRght) * log(alpha_t_right);
	sumRghtSqr += -(double)(Nplus_right) * alpha_t_right;

	if(!sumRght == 0)
		sumRghtSqr += -(double)(Nplus_right) * log(1. - exp(-(double)(sumRght) *
															const_right));

	delta = sumLeftSqr + sumRghtSqr;

	return delta;
}

double getCustomSplitNonHomogeneous(
	unsigned int n,
	char *membership,
	double *time,
	double *event,

	unsigned int eventTypeSize,
	unsigned int eventTimeSize,
	double *eventTime,

	double *response,
	double mean,
	double variance,
	unsigned int maxLevel,

	double **feature,
	unsigned int featureCount)
{
    int t;
    double delta;
    //nt = total number of events
	int *nt_left, *nt_right;
	nt_left = (int*)S_alloc(eventTypeSize, sizeof(int));
	nt_right = (int*)S_alloc(eventTypeSize, sizeof(int));

	double sum_right, sum_left;

	double *alpha_t_left, *alpha_t_right;
	alpha_t_left = (double*)S_alloc(eventTypeSize, sizeof(double));
	alpha_t_right = (double*)S_alloc(eventTypeSize, sizeof(double));

	int i, ii;
	int n_left, n_right;
	int sum_n_left, sum_n_right;
	double const_left, const_right;

	delta = 0.0;
	n_left = n_right = 0;

	sum_left = sum_right = 0;

	for (t = 0; t < eventTypeSize; t++)
		nt_left[t] = nt_right[t] = 0;


	for(i = 1; i <= n; i++)
	{
		if (membership[i] == LEFT)
		{
			nt_left[(int)event[i]] += (int)time[i];
			n_left += 1;
		}
		else
		{
			nt_right[(int)event[i]] += (int)time[i];
			n_right += 1;
		}
	}

	n_left = n_left / eventTypeSize;
	n_right = n_right / eventTypeSize;


	sum_n_left = 0;
	for(t = 0; t < eventTypeSize; t++)
		sum_n_left += nt_left[t];
	const_left = new_raph(n_left, sum_n_left);

	for(t = 0; t < eventTypeSize; t++)
		alpha_t_left[t] = (double)(nt_left[t]) * const_left;

	sum_n_right = 0;
	for(t = 0; t < eventTypeSize; t++)
		sum_n_right += nt_right[t];
	const_right = new_raph(n_right, sum_n_right);

	for(t = 0; t < eventTypeSize; t++)
		alpha_t_right[t] = (double)(nt_right[t]) * const_right;


	for(t = 0; t < eventTypeSize; t++)
		if (alpha_t_left[t] > 0.000001)
			sum_left += (double)(nt_left[t]) * log(alpha_t_left[t]);

	for (t = 0; t < eventTypeSize; t++)
		sum_left += -(double)(n_left) * alpha_t_left[t];

	if(!sum_n_left == 0)
		sum_left += -(double)(n_left) * log(1. - exp(-(double)(sum_n_left) *
													 const_left));



	for(t = 0; t < eventTypeSize; t++)
		if(alpha_t_right[t] > 0.000001)
			sum_right += (double)(nt_right[t]) * log(alpha_t_right[t]);

	for(t = 0; t < eventTypeSize; t++)
		sum_right += -(double)(n_right) * alpha_t_right[t];

	if(!sum_n_right == 0)
		sum_right += -(double)(n_right) * log(1. -
											  exp(-(double)(sum_n_right) *
												  const_right));

	delta = sum_left + sum_right;

	// Rprintf("delta=%f\n", delta);
	return delta;
}

double getCustomSplitStatisticMultivariateRegression (
	unsigned int n,
	char *membership,
	double *time,
	double *event,

	unsigned int eventTypeSize,
	unsigned int eventTimeSize,
	double *eventTime,

	double *response,
	double mean,
	double variance,
	unsigned int maxLevel,

	double **feature,
	unsigned int featureCount)
{
	// EXAMPLE:  Multivariate Regression

	// Local variables needed for this example:
	double sumLeft, sumRght;
	double sumLeftSqr, sumRghtSqr;
	double delta;

	// Left and right normalization sizes.
	unsigned int leftSize, rghtSize;

	unsigned int i;

	// Initialization of local variables:
	sumLeft = sumRght = 0.0;
	leftSize = rghtSize = 0;

	delta = 0.0;

	// In general, calculating a split statistic will require iterating
	// over all members in the parent node, and ascertaining daughter
	// membership, and performing a well defined calculation based on
	// membership.  In this example, the sum of the difference from the
	// mean for the y-outcome in each daughter node is calculated.

	for (i = 1; i <= n; i++)
	{
		// Membership will be either LEFT or RIGHT.
		if (membership[i] == LEFT)
		{
			// Add the left member to the sum.
			sumLeft += response[i] - mean;
			leftSize++;
		}
		else
		{
			// Add the right member to the sum.
			sumRght += response[i] - mean;
			rghtSize++;
		}
	}

	// Finally, we calculate the composite mean square error for each daughter.
	sumLeftSqr = pow(sumLeft, 2.0) / ((double) leftSize * variance);
	sumRghtSqr = pow(sumRght, 2.0) / ((double) rghtSize * variance);

	delta = sumLeftSqr + sumRghtSqr;

	return delta;
}

double getCustomSplitStatisticMultivariateClassification (
	unsigned int n,
	char        *membership,
	double      *time,
	double      *event,

	unsigned int eventTypeSize,
	unsigned int eventTimeSize,
	double      *eventTime,

	double      *response,
	double       mean,
	double       variance,
	unsigned int maxLevel,

	double     **feature,
	unsigned int featureCount)
{

	// EXAMPLE:  Multivariate Classification

	// Local variables needed for this example:
	unsigned int *leftClassProp, *rghtClassProp;
	double sumLeftSqr, sumRghtSqr;
	double delta;

	// Arrays needed for multivariate regression:
	double sumLeft, sumRght;

	// Left and right normalization sizes.
	unsigned int leftSize, rghtSize;

	unsigned int i, p;

	// Initialization of local variables:
	sumLeft = sumRght = 0.0;
	leftSize = rghtSize = 0;

	delta = 0.0;

	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//
	// Vector for class counts.
	leftClassProp = alloc_uivector(maxLevel);
	rghtClassProp = alloc_uivector(maxLevel);
	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//

	// Initialize the class counts.
	for (p = 1; p <= maxLevel; p++)
	{
		leftClassProp[p] = rghtClassProp[p] = 0;
	}

	// In general, calculating a split statistic will require iterating
	// over all members in the parent node, and ascertaining daughter
	// membership, and performing a well defined calculation based on
	// membership.  In this example, the left and right class counts for
	// the y-outcome in each daughter node are calculated.
	for (i = 1; i <= n; i++)
	{
		// Membership will be either LEFT or RIGHT.  Package specific
		// constants are contained in global.h.
		if (membership[i] == LEFT)
		{
			// Add the left member to the left class proportion.
			leftClassProp[(unsigned int) response[i]]++;
			leftSize++;
		}
		else
		{
			// Add the right member to the left class proportion.
			rghtClassProp[(unsigned int) response[i]]++;
			rghtSize++;
		}
	}

	for (p = 1; p <= maxLevel; p++)
	{
		sumLeft += pow((double) leftClassProp[p], 2.0);
		sumRght += pow((double) rghtClassProp[p], 2.0);
	}

	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//
	// Vector for class counts.
	dealloc_uivector(leftClassProp, maxLevel);
	dealloc_uivector(rghtClassProp, maxLevel);
	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//

	sumLeftSqr = sumLeft / leftSize;
	sumRghtSqr = sumRght / rghtSize;

	delta = sumLeftSqr + sumRghtSqr;

	return delta;
}

double getCustomSplitStatisticSurvival (
	unsigned int n,
	char        *membership,
	double      *time,
	double      *event,

	unsigned int eventTypeSize,
	unsigned int eventTimeSize,
	double      *eventTime,

	double      *response,
	double       mean,
	double       variance,
	unsigned int maxLevel,

	double     **feature,
	unsigned int featureCount)
{
	// EXAMPLE:  Survival (logrank)

	// Local variables needed for this example:
	unsigned int *nodeLeftEvent, *nodeParentEvent;
	unsigned int *nodeLeftAtRisk, *nodeParentAtRisk;

	double delta, deltaNum, deltaDen;
	unsigned int i, k;

	// Initialization of local variables:
	deltaNum = deltaDen = 0.0;

	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//
	// Vector for event counts.
	nodeLeftEvent = alloc_uivector(eventTimeSize);
	nodeParentEvent = alloc_uivector(eventTimeSize);

	// Vector for at risk counts.
	nodeLeftAtRisk = alloc_uivector(eventTimeSize);
	nodeParentAtRisk  = alloc_uivector(eventTimeSize);

	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//

	// Reset the node specific counts needed for the statistic.
	for (k = 1; k <= eventTimeSize; k++)
	{
		nodeParentEvent[k] = 0;
		nodeParentAtRisk[k] = 0;
		nodeLeftEvent[k] = 0;
		nodeLeftAtRisk[k] = 0;
	}

	k = eventTimeSize;
	i = n;

	// Iterate over all individuals.  Note that they arrive sorted by time
	// in increasing order.  We parse them in decreasing order.
	while ((i > 0) && (k > 0))
	{
		// Initialize the parent event count.
		if (eventTime[k] <= time[i])
		{
			// The member is still at risk!
			nodeParentAtRisk[k]++;

			// Membership will be either LEFT or RIGHT.
			if (membership[i] == LEFT)
			{
				nodeLeftAtRisk[k]++;
			}

			// Did the member experience an event?
			if (eventTime[k] == time[i])
			{
				if (event[i] > 0)
				{
					nodeParentEvent[k]++;

					// Membership will be either LEFT or RIGHT.
					if (membership[i] == LEFT)
					{
						nodeLeftEvent[k]++;
					}
				}
			}

			// Examine the previous individual.
			i--;
		}
		else
		{
			// Examine the previous event time.
			k--;
		}
	}

	// Adjust the at risk counts to achieve the step function.
	for (k = eventTimeSize; k > 1; k--)
	{
		nodeParentAtRisk[k - 1] = nodeParentAtRisk[k] + nodeParentAtRisk[k - 1];
		nodeLeftAtRisk[k - 1] = nodeLeftAtRisk[k] + nodeLeftAtRisk[k - 1];
	}

	// Iterate over the distinct event times and acquire the numerator and
	// denominator of the test.
	for (k = 1; k <= eventTimeSize; k++) {
		deltaNum = deltaNum + ((double) nodeLeftEvent[k] -
							   ((double) (nodeLeftAtRisk[k] *
							   			  nodeParentEvent[k]) /
							  	nodeParentAtRisk[k]));

		// Log-Rank denominator requires that there be at least two at risk.
		if (nodeParentAtRisk[k] >= 2)
		{
			deltaDen = deltaDen + (((double) nodeLeftAtRisk[k] /
									nodeParentAtRisk[k]) *
								  (1.0 - ((double) nodeLeftAtRisk[k] /
								   nodeParentAtRisk[k])) *
								  ((double) (nodeParentAtRisk[k] -
								   			 nodeParentEvent[k]) /
								   (nodeParentAtRisk[k] - 1)) *
								  nodeParentEvent[k]);
		}
	}


	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//
	// Vector for event counts.
	dealloc_uivector(nodeLeftEvent, eventTimeSize);
	dealloc_uivector(nodeParentEvent, eventTimeSize);

	// Vector for at risk counts.
	dealloc_uivector(nodeLeftAtRisk, eventTimeSize);
	dealloc_uivector(nodeParentAtRisk, eventTimeSize);
	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//

	deltaNum = fabs(deltaNum);
	deltaDen = sqrt(deltaDen);
	if (deltaDen <= 1.0e-9)
	{
		if (deltaNum <= 1.0e-9)
		{
			delta = 0.0;
		}
		else
		{
			delta = deltaNum / deltaDen;
		}
	}
	else
	{
		delta = deltaNum / deltaDen;
	}

	return delta;
}

double getCustomSplitStatisticCompetingRisk (
	unsigned int n,
	char        *membership,
	double      *time,
	double      *event,

	unsigned int eventTypeSize,
	unsigned int eventTimeSize,
	double      *eventTime,

	double      *response,
	double       mean,
	double       variance,
	unsigned int maxLevel,

	double     **feature,
	unsigned int featureCount)
{
	// EXAMPLE:  Competing Risk (logrankCR)

	// Local variables needed for this example:
	unsigned int *nodeLeftEvent, *nodeParentEvent;
	unsigned int *nodeLeftAtRisk, *nodeParentAtRisk;

	// Local CR variable needed for this example:
	unsigned int **nodeLeftEventCR, **nodeParentEventCR;
	unsigned int **nodeLeftInclusiveAtRisk, **nodeParentInclusiveAtRisk;

	double delta, deltaNum, deltaSubNum, deltaDen, deltaSubDen;
	unsigned int i, j, k, r, s;

	// Initialization of local variables:
	deltaNum = deltaDen = 0.0;

	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//
	// Vector for event counts.
	nodeLeftEvent = alloc_uivector(eventTimeSize);
	nodeParentEvent = alloc_uivector(eventTimeSize);

	// Vector for at risk counts.
	nodeLeftAtRisk = alloc_uivector(eventTimeSize);
	nodeParentAtRisk = alloc_uivector(eventTimeSize);
	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//

	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//
	// Matrix containing event counts at each event time.
	nodeParentEventCR = alloc_uimatrix(eventTypeSize, eventTimeSize);
	nodeLeftEventCR = alloc_uimatrix(eventTypeSize, eventTimeSize);

	// Matrix containing event inclusive at risk counts at each event time.
	nodeParentInclusiveAtRisk = alloc_uimatrix(eventTypeSize, eventTimeSize);
	nodeLeftInclusiveAtRisk = alloc_uimatrix(eventTypeSize, eventTimeSize);
	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//

	// Reset the node specific counts needed for the statistic.
	for (k = 1; k <= eventTimeSize; k++)
	{
		nodeParentEvent[k] = 0;
		nodeParentAtRisk[k] = 0;
		nodeLeftEvent[k] = 0;
		nodeLeftAtRisk[k] = 0;

		for (j = 1; j <= eventTypeSize; j++)
		{
			nodeParentEventCR[j][k] = 0;
			nodeLeftEventCR[j][k] = 0;
			nodeParentInclusiveAtRisk[j][k] = 0;
			nodeLeftInclusiveAtRisk[j][k] = 0;
		}
	}

	k = eventTimeSize;
	i = n;

	// Iterate over all individuals.  Note that they arrive sorted by time
	// in increasing order.  We parse them in decreasing order.
	while ((i > 0) && (k > 0))
	{
		// Initialize the parent event count.
		if (eventTime[k] <= time[i])
		{
			// The member is still at risk!
			nodeParentAtRisk[k]++;

			// Membership will be either LEFT or RIGHT.
			if (membership[i] == LEFT)
			{
				nodeLeftAtRisk[k]++;
			}

			// Did the member experience an event?
			if (eventTime[k] == time[i])
			{
				if (event[i] > 0)
				{
					nodeParentEventCR[(unsigned int) event[i]][k]++;
					nodeParentEvent[k]++;

					// Membership will be either LEFT or RIGHT.
					if (membership[i] == LEFT)
					{
						nodeLeftEventCR[(unsigned int) event[i]][k]++;
					}
				}
			}

			// Examine the previous individual.
			i--;
		}
		else
		{
			// Examine the previous event time.
			k--;
		}
	}

	// Adjust the at risk counts to achieve the step function.
	for (k = eventTimeSize; k > 1; k--)
	{
		nodeParentAtRisk[k - 1] = nodeParentAtRisk[k] + nodeParentAtRisk[k - 1];
		nodeLeftAtRisk[k - 1] = nodeLeftAtRisk[k] + nodeLeftAtRisk[k - 1];
	}

	// Finalize the left and right inclusive at risk counts.
	for (k = 1; k <= eventTimeSize; k++)
	{
		for (j = 1; j <= eventTypeSize; j++)
		{
			nodeParentInclusiveAtRisk[j][k] = nodeParentAtRisk[k];
			nodeLeftInclusiveAtRisk[j][k] = nodeLeftAtRisk[k];
			for (s = 1; s < k; s++)
			{
				for (r = 1; r <= eventTypeSize; r++)
				{
					if (j != r)
					{
						nodeParentInclusiveAtRisk[j][k] += nodeParentEventCR[r][s];
						nodeLeftInclusiveAtRisk[j][k] += nodeLeftEventCR[r][s];
					}
				}
			}
		}
	}

	// Iterate over the event types and distinct event times and acquire the
	// numerator and denominator of the test.
	for (j = 1; j <= eventTypeSize; j++)
	{
		deltaSubNum = 0;
		for (k = 1; k <= eventTimeSize; k++)
		{
			deltaSubNum = deltaSubNum +
						  (nodeLeftEventCR[j][k] -
						  (nodeParentEventCR[j][k] *
						   ((double) nodeLeftInclusiveAtRisk[j][k] /
						   	nodeParentInclusiveAtRisk[j][k])));
		}
		deltaNum = deltaNum + deltaSubNum;

		deltaSubDen = 0;
		for (k = 1; k <= eventTimeSize; k++)
		{
			// Log-Rank CR denominator requires that there be at least two at risk.
			if (nodeParentAtRisk[k] >= 2)
			{
				deltaSubDen = deltaSubDen +
							  ((nodeParentEventCR[j][k] *
							  ((double) nodeLeftInclusiveAtRisk[j][k] /
							   nodeParentInclusiveAtRisk[j][k])) *
							  (1.0 - ((double) nodeLeftInclusiveAtRisk[j][k] /
							   nodeParentInclusiveAtRisk[j][k])) *
							  ((double) (nodeParentInclusiveAtRisk[j][k] -
							   nodeParentEventCR[j][k]) /
							  (nodeParentInclusiveAtRisk[j][k] - 1)));
			}
		}
		deltaDen = deltaDen + deltaSubDen;
	}

	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//
	// Vector for event counts.
	dealloc_uivector(nodeLeftEvent, eventTimeSize);
	dealloc_uivector(nodeParentEvent, eventTimeSize);

	// Vector for at risk counts.
	dealloc_uivector(nodeLeftAtRisk, eventTimeSize);
	dealloc_uivector(nodeParentAtRisk, eventTimeSize);
	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//

	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//
	// Matrix containing event counts at each event time.
	dealloc_uimatrix(nodeParentEventCR, eventTypeSize, eventTimeSize);
	dealloc_uimatrix(nodeLeftEventCR, eventTypeSize, eventTimeSize);

	// Matrix containing event inclusive at risk counts at each event time.
	dealloc_uimatrix(nodeParentInclusiveAtRisk, eventTypeSize, eventTimeSize);
	dealloc_uimatrix(nodeLeftInclusiveAtRisk, eventTypeSize, eventTimeSize);
	//---------- *** WARNING MEMORY MANIPULATION WARNING *** -----------//

	deltaNum = fabs(deltaNum);
	deltaDen = sqrt(deltaDen);
	if (deltaDen <= 1.0e-9)
	{
		if (deltaNum <= 1.0e-9)
		{
			delta = 0.0;
		}
		else
		{
			delta = deltaNum / deltaDen;
		}
	}
	else
	{
		delta = deltaNum / deltaDen;
	}

	return delta;
}

/*
  Memory allocation and deallocation.
    [nh] = the length of the array

  Note that indexing is one-based:
    array[1] ... array[nh]

  Multi-dimensional array allocationis
  accomplished via multiple one-dimensional
  array allocations.
*/
unsigned int *alloc_uivector(unsigned int nh)
{
	return (unsigned int*) malloc((size_t) ((nh + 1) * (sizeof(unsigned int))));
}

void dealloc_uivector(unsigned int *v, unsigned int nh)
{
	free((char *) v);
}

double *alloc_dvector(double *v, unsigned int nh)
{
	return (double *) malloc((size_t) ((nh + 1) * (sizeof(double))));
}

void dealloc_dvector(double *v, unsigned int nh)
{
	free((char *) v);
}

unsigned int **alloc_uimatrix(unsigned int n2h, unsigned int nh)
{
	unsigned int **v = (unsigned int **) malloc((size_t) ((n2h + 1) *
					   (sizeof(unsigned int *))));

	for (unsigned int i = 1; i <= n2h; i++)
	{
		v[i] = alloc_uivector(nh);
	}

	return v;
}

void dealloc_uimatrix(unsigned int **v, unsigned int n2h, unsigned int nh)
{
	for (unsigned int i = 1; i <= n2h; i++)
	{
		dealloc_uivector(v[i], nh);
	}

	free((char *) v);
}