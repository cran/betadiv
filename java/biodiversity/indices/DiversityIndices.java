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

import java.util.HashMap;
import java.util.Map;

public class DiversityIndices {

	private final Map<DiversityIndex, Object> innerMap;
	
	public DiversityIndices() {
		super();
		innerMap = new HashMap<DiversityIndex, Object>();
	}

	protected static enum DiversityIndex {Alpha,
		Gamma,
		Beta
	}
	
	public static enum BetaIndex {
		Simpson, 
		Sorensen,
		Nestedness
	}  
	
	protected void setAlphaDiversity(double alpha) {
		innerMap.put(DiversityIndex.Alpha, alpha);
	}
	
	public double getAlphaDiversity() {
		if (innerMap.containsKey(DiversityIndex.Alpha)) {
			return (Double) innerMap.get(DiversityIndex.Alpha);
		} else {
			return -1d;
		}
	}
	
	protected void setGammaDiversity(double gamma) {
		innerMap.put(DiversityIndex.Gamma, gamma);
	}

	public double getGammaDiversity() {
		if (innerMap.containsKey(DiversityIndex.Gamma)) {
			return (Double) innerMap.get(DiversityIndex.Gamma);
		} else {
			return -1d;
		}
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private Map<BetaIndex, Double> getBetaIndices() {
		if (!innerMap.containsKey(DiversityIndex.Beta)) {
			innerMap.put(DiversityIndex.Beta, new HashMap<BetaIndex, Double>());
		}
		return (Map) innerMap.get(DiversityIndex.Beta);
	}
	
	protected void setBetaDiversity(BetaIndex ind, double value) {
		getBetaIndices().put(ind, value);
	}

	public double getBetaDiversity(BetaIndex ind) {
		if (getBetaIndices().containsKey(ind)) {
			return getBetaIndices().get(ind);
		} else {
			return -1d;
		}
	}
	
}
