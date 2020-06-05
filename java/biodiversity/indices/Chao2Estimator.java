/*
 * This file is part of the betadiversityindices library
 *
 * Author Mathieu Fortin - Canadian Forest Service
 * Copyright (C) 2020 Her Majesty the Queen in right of Canada
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */package biodiversity.indices;

import java.security.InvalidParameterException;
import java.util.List;
import java.util.Map;

import biodiversity.indices.IndexUtility.ValidatedHashMap;
import biodiversity.indices.MultipleSiteIndex.SpeciesFreqMap;
import repicea.math.Matrix;
import repicea.stats.estimates.SimpleEstimate;

/**
 * This class computes the estimator of species richness developed known as Chao2. 
 *  
 * @author Mathieu Fortin
 * @see Chao, A., & Colwell, R. K. (2017). Thirty years of progeny from Chao’s inequality: Estimating and comparing richness with incidence data
 * and incomplete sampling. SORT: Statistics and Operations Research Transactions, 41(1), 3–54.
 */
public class Chao2Estimator {

	
	/**
	 * Returns the estimate of species richness in the context of sampling without replacement. 
	 * For sampling with replacement, use the getChao2Estimate(Map<String, List> sample) method instead.
	 * If the populationSize parameter is set to 0, the method works under the assumption of
	 * random sampling with replacement.
	 * @param sample a Map instance that stands for the sample
	 * @param populationSize an integer >= 0. 
	 * @return a SimpleEstimate instance
	 */
	@SuppressWarnings("rawtypes")
	public static SimpleEstimate getChao2Estimate(Map<String, List> sample, int populationSize) {
		if (populationSize < 0) {
			throw new InvalidParameterException("The population size must be equal to or greater than 0!");
		}
		if (sample == null || sample.size() < 2) {
			throw new InvalidParameterException("The sample parameter must have at least two plots or sites!");
		}
		ValidatedHashMap<String, List> validatedMap = IndexUtility.validateMap(sample);
		boolean withReplacement = populationSize == 0;
		
		int nbPlots = validatedMap.size();

		SpeciesFreqMap speciesFreqMap = getSpeciesFreqMap(validatedMap);
		int f1 = speciesFreqMap.getNbSpeciesWithThisFreq(1);
		int f2 = speciesFreqMap.getNbSpeciesWithThisFreq(2);
		int s = IndexUtility.getUniqueSpeciesList(sample).size();
		double f1_f2 = ((double) f1) / f2;
		double k = ((double) (nbPlots - 1)) / nbPlots;
		double chao2;
		double variance;
		if (f2 == 0) {
			chao2 = s + ((double) (nbPlots - 1)) / nbPlots * f1 * (f1 - 1) / 2; 
			variance = k * f1 * (f1 - 1) / 2d + k * k * f1 * (2 * f1 - 1) * (2 * f1 - 1) / 4d + k * k * f1 * f1 * f1 * f1  / (4d * chao2); 
		} else {
			double w = ((double) nbPlots) / (nbPlots - 1); 
			if (withReplacement) {
				chao2 = s +  f1 * f1 / (w * 2d * f2); 
				variance = f2 * (.5 * k * f1_f2 * f1_f2 + k * k * f1_f2 * f1_f2 * f1_f2 + .25 * k * k * f1_f2 * f1_f2 * f1_f2 * f1_f2); 
			} else {
				double r = ((double) nbPlots) / populationSize;
				chao2 = s +  f1 * f1 / (w * 2d * f2 + r / (1 - r) * f1); 
				double f0 = chao2 - s;
				double f0_f1 = f0 / f1;
				double varianceInnerTerm = 2 * w * f2 * f0 * f0 + f1 * f1 * f0;
				variance = f0 + varianceInnerTerm * varianceInnerTerm / (f1 * f1 * f1 * f1 * f1) +	4 * w * w * f2 * (f0_f1 * f0_f1 * f0_f1 * f0_f1);
			}
		}
		SimpleEstimate estimate = new SimpleEstimate();
		Matrix mean = new Matrix(1,1);
		mean.m_afData[0][0] = chao2;
		Matrix var = new Matrix(1,1);
		var.m_afData[0][0] = variance;
		estimate.setMean(mean);
		estimate.setVariance(var);
		return estimate;
	}

	
	/**
	 * Returns the estimate of species richness in the context of sampling with replacement. 
	 * For sampling without replacement, use the getChao2Estimate(Map<String, List> sample,
	 * int populationSize) method instead.
	 * @param sample a Map instance that stands for the sample
	 * @return a SimpleEstimate instance
	 */
	public static SimpleEstimate getChao2Estimate(Map<String, List> sample) {
		return getChao2Estimate(sample, 0);
	}
	
	@SuppressWarnings("rawtypes")
	private static SpeciesFreqMap getSpeciesFreqMap(ValidatedHashMap<String, List> vMap) {
		SpeciesFreqMap speciesMap = new SpeciesFreqMap();
		for (List speciesList : vMap.values()) {
			for (Object s : speciesList) {
				if (!speciesMap.containsKey(s)) {
					speciesMap.put(s, 0);
				}
				speciesMap.put(s, speciesMap.get(s) + 1);
			}
		}
		return speciesMap;
	}

}
