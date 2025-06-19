package org.example.GA;

public class Configuration {
    private final String instance;
    private final int populationSize;
    private final String selectionMethod;
    private final double elitismRate;
    private final int TSRate;
    private final String mutationMethod;
    private final double mutationRate;
    private final String crossoverMethod;
    private final double crossoverRate;
    private final int numberOfElites;
    private final int LSRate;

    public Configuration(String instance, int populationSize, String selectionMethod, double elitismRate, int TSRate, String mutationMethod, double mutationRate, String crossoverMethod, double crossoverRate, int numberOfElites, int LSRate) {
        this.instance = instance;
        this.populationSize = populationSize;
        this.selectionMethod = selectionMethod;
        this.elitismRate = elitismRate;
        this.TSRate = TSRate;
        this.mutationMethod = mutationMethod;
        this.mutationRate = mutationRate;
        this.crossoverMethod = crossoverMethod;
        this.crossoverRate = crossoverRate;
        this.numberOfElites = numberOfElites;
        this.LSRate = LSRate;
    }
    public String getInstance() {
        return instance;
    }

    public int getPopulationSize() {
        return populationSize;
    }

    public String getSelectionMethod() {
        return selectionMethod;
    }

    public double getElitismRate() {
        return elitismRate;
    }

    public int getTSRate() {
        return TSRate;
    }

    public String getMutationMethod() {
        return mutationMethod;
    }

    public double getMutRate() {
        return mutationRate;
    }

    public String getCrossoverMethod() {
        return crossoverMethod;
    }

    public double getCrossRate() {
        return crossoverRate;
    }

    public int getNumberOfElites() {
        return numberOfElites;
    }

    public int getLSRate() {
        return LSRate;
    }
}
