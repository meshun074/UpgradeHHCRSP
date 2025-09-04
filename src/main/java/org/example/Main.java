package org.example;

import org.example.Data.InstancesClass;
import org.example.Data.ReadData;
import org.example.GA.*;

import java.io.File;
import java.io.PrintStream;


//TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
public class Main {
    public static InstancesClass instance;
    /*
    Usage: main arg1-int:parameter settings (0-65) arg2-int:problem size (0-6) arg3-int:instance number (1-10) arg4-int:seed
     */

    @SuppressWarnings("CallToPrintStackTrace")
    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Usage: java GeneticAlgorithmRunner <config-file>");
            return;
        }
        for (int i = 0; i < 1; i++) {
            long randomSeed;
            long instanceNumber;
            String instanceName;

            try {
                // Read configuration file
                File configFile = new File(args[0]);
                if (args[0].contains("HHCRSP")) {
                    Config1 config1 = Config1.read(configFile);
                    instanceNumber = config1.getInstanceIndex();
                    instanceName = config1.getInstanceName();
                    int runCount = Integer.parseInt(args[1]);
                    randomSeed = System.currentTimeMillis() + runCount;
                    //create result directory
                    String resultDir = "src/main/java/org/example/Config_" + instanceName + "_" + instanceNumber + "_results";
                    new File(resultDir).mkdirs();
                    // Read dataset
                    PrintStream fileout = new PrintStream(resultDir + "/Result_" + instanceName + "_" + instanceNumber + "_" + runCount + "_" + "_" + randomSeed + ".txt");
                    System.setOut(fileout);
                    System.out.printf("Config : instanceNumber=%d, seed=%d\n", instanceNumber, randomSeed);

                    instance = ReadData.read(new File("src/main/java/org/example/Data/kummer/" + instanceName));

                } else {
                    Config config = Config.read(configFile);

                    // Extract parameters from JSON
                    int problemSize = config.getProblemSize();
                    instanceNumber = config.getInstanceIndex();

                    int runCount = Integer.parseInt(args[1]);
                    randomSeed = System.currentTimeMillis() + runCount;

                    String[] Instances = {"10", "25", "50", "75", "100", "200", "300"};
                    instanceName = Instances[problemSize];

                    //create result directory
                    String resultDir = "src/main/java/org/example/BSTest_" + problemSize + "_" + instanceNumber + "_results";
                    new File(resultDir).mkdirs();

                    // Read dataset
                    PrintStream fileout = new PrintStream(resultDir + "/Result_" + instanceName + "_" + instanceNumber + "_" + runCount + "_" + randomSeed + ".txt");
                    System.setOut(fileout);
                    System.out.printf("Config Parameters: ProblemSize=%d, instanceNumber=%d, seed=%d\n", problemSize, instanceNumber, randomSeed);

                    instance = ReadData.read(new File("src/main/java/org/example/Data/instance/" + instanceName + "_" + instanceNumber + ".json"));
                }

                Configuration config = parseArguments.getConfiguration(args);

                GeneticAlgorithm ga = new GeneticAlgorithm(config, 1000, instance);
                Chromosome bestChromosome = ga.start();

                assert bestChromosome != null;
                System.out.println("----------------- Solution ----------------------");
                System.out.println("Instance_" + instanceName + "_" + instanceNumber + " Best Fitness: " + Math.round(bestChromosome.getFitness() * 1000.0) / 1000.0);
                System.out.println("Total Distance: " + bestChromosome.getTotalTravelCost() + " Total Tardiness: " + bestChromosome.getTotalTardiness() + " Highest Tardiness: " + bestChromosome.getHighestTardiness());
                bestChromosome.showSolution(0);
            } catch (Exception e) {
                e.printStackTrace();
            }


        }

    }
}