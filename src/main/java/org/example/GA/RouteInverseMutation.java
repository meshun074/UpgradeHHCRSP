package org.example.GA;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class RouteInverseMutation implements Runnable {
    @Override
    public void run() {
        ga.getMutationChromosomes().add(mutation());
    }

    private final GeneticAlgorithm ga;
    private final Random rand;
    private final Chromosome ch;

    public RouteInverseMutation(GeneticAlgorithm ga, Chromosome ch) {
        this.ga = ga;
        this.rand = new Random(System.nanoTime());
        this.ch = ch;
    }

    public Chromosome mutation() {
        Chromosome mutant = new Chromosome(ch.getGenes(), ch.getFitness(), true);
        List<Integer>[] genes = mutant.getGenes();
        int selectedRoute =0;
        List<Integer> route = null;
        int end =3;
        do{
            end--;
            if(route!=null){
                genes[selectedRoute] = route;
            }
            //Select random route to mutate
            selectedRoute = rand.nextInt(genes.length);
            route = genes[selectedRoute];
            int routeSize = route.size();

            // Only mutate if route has at least 2 elements
            if(routeSize > 1){
                List<Integer> newRoute = new ArrayList<>(route);
                // Apply mutation based on route size
                if (routeSize == 2) {
                    // Simple swap for size 2
                    Collections.reverse(newRoute);
                }else {
                    int index = rand.nextInt(routeSize);

                    if (routeSize - index == 1) {
                        // Handle end of route cases
                        swap(newRoute, index, rand.nextBoolean() ? index - 1 : index - 2);
                    } else if (routeSize - index == 2) {
                        swap(newRoute, index, index + 1);
                    } else {
                        swap(newRoute, index, rand.nextBoolean() ? index + 1 : index + 2);
                    }
                }

                //Update the route
                genes[selectedRoute] = newRoute;
            }

            EvaluationFunction.Evaluate(mutant);
        }while (mutant.getFitness() == Double.POSITIVE_INFINITY&& end>0);
        return mutant;
    }
    private void swap(List<Integer> newRoute, int index1, int index2) {
        int key1 = newRoute.get(index1);
        newRoute.set(index1, newRoute.get(index2));
        newRoute.set(index2, key1);
    }
}
