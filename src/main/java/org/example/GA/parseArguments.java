package org.example.GA;

public class parseArguments {
    public static Configuration getConfiguration(String[] args) {
        String instance = "";
        int populationSize = 0;
        String selectionMethod = "";
        double elitismRate = 0;
        int TSRate = 0;
        String mutationMethod = "";
        double mutRate = 0;
        String crossoverMethod = "";
        double crossRate =0;
        int numberOfElites= 0;
        int LSRate =0;
        boolean LSStart = true;
        double LSStartRate = 0.0;
//        String arg = "--popSize 300 --mutRate 0.05 --mutMethod S --crossRate 1.0 --crossMethod BS --selection RW --elitism 0.1 --numberOfElites 5 --LSRate 10 --TSRate 4 --instance src/main/java/org/example/Data/instance/200_1.json";
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--instance":
                    instance = args[++i];
                    break;
                case "--popSize":
                    populationSize = Integer.parseInt(args[++i]);
                    break;
                case "--mutRate":
                    mutRate = Double.parseDouble(args[++i]);
                    break;
                case "--mutMethod":
                    mutationMethod = args[++i];
                    break;
                case "--crossRate":
                    crossRate = Double.parseDouble(args[++i]);
                    break;
                case "--crossMethod":
                    crossoverMethod = args[++i];
                    break;
                case "--selection":
                    selectionMethod = args[++i];
                    break;
                case "--elitism":
                    elitismRate = Double.parseDouble(args[++i]);
                    break;
                case "--numberOfElites":
                    numberOfElites = Integer.parseInt(args[++i]);
                    break;
                case "--LSRate":
                    LSRate = Integer.parseInt(args[++i]);
                    break;
                case "--LSStart":
                    LSStart = Boolean.parseBoolean(args[++i]);
                    break;
                case "--LSStartRate":
                    LSStartRate = Double.parseDouble(args[++i]);
                    break;
                case "--TSRate":
                    TSRate = Integer.parseInt(args[++i]);
                    break;
            }
        }
        return new Configuration(instance,populationSize,selectionMethod,elitismRate,TSRate,mutationMethod,mutRate,crossoverMethod,crossRate,numberOfElites,LSRate,LSStart,LSStartRate);
    }
}
