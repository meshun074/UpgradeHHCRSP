package org.example.GA;

import org.example.Data.Caregiver;
import org.example.Data.InstancesClass;
import org.example.Data.Patient;
import org.example.Data.Required_Caregiver;

import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

public class EvaluationFunction {
    private static InstancesClass dataset;
    private static double[][] distances;
    private static Patient[] allPatients;
    private static Caregiver[] allCaregivers;
    private static final ThreadLocal<Set<Integer>> trackHolder = ThreadLocal.withInitial(() -> new HashSet<>(150));


    public static void initialize(InstancesClass data) {
        dataset = data;
        distances = dataset.getDistances();
        allPatients = dataset.getPatients();
        allCaregivers = dataset.getCaregivers();
    }

    public static void EvaluateFitness(List<Chromosome> population) {

                try (ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors())) {
                    try {
                        List<Future<?>> futures = population.stream()
                                .map(ch -> executor.submit(() -> Evaluate(ch)))
                                .collect(Collectors.toList());
                        for (Future<?> f : futures) f.get(); // Wait for completion
                    } catch (ExecutionException | InterruptedException e) {
                        throw new RuntimeException(e);
                    } finally {
                        executor.shutdown();
                    }
                }
    }
    private static void Evaluate(Chromosome temp, Chromosome base, boolean isInvalid){
        int[] routeEndPoint = new int[allCaregivers.length];
        Arrays.fill(routeEndPoint, -1);
        if(isInvalid){
            Evaluate(temp);
        }
        int first = temp.getFirst();
        int second = temp.getSecond();
        int firstPosition = temp.getFirstPosition();
        int secondPosition = temp.getSecondPosition();
        int size = 1;
        routeEndPoint[first] = firstPosition;
        if(second !=-1){
            size++;
            routeEndPoint[second] = secondPosition;
        }
        int[] routeMove = new int[size];
        int[] positionMove = new int[size];
        routeMove[0] = first;
        positionMove[0] = firstPosition;
        if(size > 1){
            routeMove[1] = second;
            positionMove[1] = secondPosition;
        }

        removeAffectedPatient(routeMove,positionMove, base,routeEndPoint);
        Shift[] shifts = base.getCaregiversRouteUp();
        Shift[] tempShifts = temp.getCaregiversRouteUp();

        for(int i = 0; i < routeEndPoint.length; i++){
            List<Integer> route = new ArrayList<>(shifts[i].getRoute());
            List<Double> currentTime  = new ArrayList<>(shifts[i].getCurrentTime());
            List<Double> travelCost = new ArrayList<>(shifts[i].getTravelCost());
            travelCost.remove(travelCost.size()-1);
            List<Double> tardiness = new ArrayList<>(shifts[i].getTardiness());
            List<Double> maxTardiness = new ArrayList<>(shifts[i].getMaxTardiness());
            if(routeEndPoint[i] != -1){
                int index = routeEndPoint[i];
                route.subList(index, route.size()).clear();
                index++;
                travelCost.subList(index, travelCost.size()).clear();
                currentTime.subList(index, currentTime.size()).clear();
                tardiness.subList(index, tardiness.size()).clear();
                maxTardiness.subList(index, maxTardiness.size()).clear();
                tempShifts[i] = new Shift(shifts[i].getCaregiver(),route,currentTime,travelCost,tardiness,maxTardiness);
            }else {
                tempShifts[i] = new Shift(shifts[i].getCaregiver(),route,currentTime,travelCost,tardiness,maxTardiness);
            }
        }
        double totalTravelCost = 0;
        double totalTardiness = 0;
        double highestTardiness = 0;
        for(Shift s : tempShifts){
            if( s.getTravelCost().isEmpty()){
                s.addTravelCost(0.0);
                s.addTardiness(0.0);
                s.initializeMaxTardiness(0.0);
                s.addCurrentTime(0.0);
            }
            totalTravelCost+=  s.getTravelCost().get(s.getTravelCost().size()-1);
            totalTardiness+= s.getTardiness().isEmpty()?0: s.getTardiness().get(s.getTardiness().size()-1);
            double maxTardiness = s.getMaxTardiness().isEmpty()?0: s.getMaxTardiness().get(s.getMaxTardiness().size()-1);
            highestTardiness = Math.max(highestTardiness,maxTardiness);
        }
        temp.setTotalTravelCost(totalTravelCost);
        temp.setTotalTardiness(totalTardiness);
        temp.setHighestTardiness(highestTardiness);
        temp.setFitness(0.0);
        evaluate(temp,routeEndPoint);

    }
    public static void removeAffectedPatient(int[] routeMove, int[] positionMove, Chromosome base, int[] routeEndPoint) {
        Map<Integer, Set<Integer>> patientToRoutesMap = base.getPatientToRoutesMap();
        List<Integer>[] genes = base.getGenes();
        for(int j = 0; j < routeMove.length; j++){
            int first = routeMove[j];
            int firstPosition = positionMove[j];
            List<Integer> currentRoute = genes[first];
            for (int i = firstPosition; i < currentRoute.size(); i++) {
                int patient = currentRoute.get(i);
                Patient p = allPatients[patient];
                if (p.getRequired_caregivers().length > 1) {
                    int routeIndex = getRouteIndexMethod(first, patientToRoutesMap.get(patient));
                    int patientIndex = genes[routeIndex].indexOf(patient);
                    if (routeEndPoint[routeIndex] == -1 || routeEndPoint[routeIndex] > patientIndex) {
                        routeEndPoint[routeIndex] = patientIndex;
                        int[] newRouteMove = {routeIndex};
                        int[] newPositionMove = {patientIndex};
                        removeAffectedPatient(newRouteMove,newPositionMove, base, routeEndPoint);
                    }
                }
            }
        }
    }

    private static int getRouteIndexMethod(int first, Set<Integer> routes) {
        if(routes==null) return -1; // Patient not found
        for(int route : routes){
            if(route != first){
                return route; // Return the first alternative route
            }
        }
        return -1; // other route not found
    }

    private static void evaluate(Chromosome temp, int[] routeEndPoint) {
        Shift[] shifts = temp.getCaregiversRouteUp();
        Set<Integer> track = new HashSet<>(100);
        List<Integer>[] genes = temp.getGenes();
        for(int i = 0; i < routeEndPoint.length; i++){
            List<Integer> route = genes[i];
            Shift caregiver = shifts[i];
            int routeEnd = routeEndPoint[i];
            if(routeEnd != -1){
                for(int j = routeEnd; j < route.size(); j++){
                    int patient = route.get(j);
                    if(!caregiver.getRoute().contains(patient)){
                        if(patientAssignment(temp,patient,caregiver,shifts,i,track)){
                            temp.setFitness(Double.POSITIVE_INFINITY);
                            return;
                        }
                        UpdateCost(temp);
                        track.clear();
                    }
                }
            }
        }
        for(Shift s : shifts){
            int lastLocationId = s.getRoute().isEmpty()? 0:s.getRoute().get(s.getRoute().size() - 1) + 1;
            double travelCost = distances[lastLocationId][0];
            temp.updateTotalTravelCost(travelCost);
            s.updateTravelCost(travelCost);
        }
        UpdateCost(temp);
    }

    public static void Evaluate(Chromosome ch) {
        Shift[] routes = ch.getCaregiversRouteUp();
        for (int i = 0; i < routes.length; i++) {
            routes[i].resetShift();
        }
        ch.setCaregiversRouteUp(routes);
        ch.setHighestTardiness(0);
        ch.setTotalTardiness(0);
        ch.setTotalTravelCost(0);
        ch.setFitness(Double.POSITIVE_INFINITY);

        Set<Integer> track = trackHolder.get(); // Initial capacity for small patient sets
        List<Integer>[] genes = ch.getGenes();

        try {
            for (int i = 0; i < routes.length; i++) {
                List<Integer> route = genes[i];
                Shift caregiver1 = routes[i];

                for (int patient : route) {
                    if (!caregiver1.getRoute().contains(patient)) {
                        if (patientAssignment(ch, patient, caregiver1, routes, i, track)) {
                            ch.setFitness(Double.POSITIVE_INFINITY);
                            return;
                        }
                        track.clear();
                    }
                }
            }
        }
        finally {
            track.clear();
        }

        // Final cost calculations
        for (Shift s : routes) {
            int lastLocationId = s.getRoute().isEmpty()? 0:s.getRoute().get(s.getRoute().size() - 1) + 1;
            double returnCost = distances[lastLocationId][0];
            ch.updateTotalTravelCost(returnCost);
            s.updateTravelCost(returnCost);
        }

        UpdateCost(ch);
    }

    public static void EvaluateNew(Chromosome ch) {
        ch.setHighestTardiness(0);
        ch.setTotalTardiness(0);
        ch.setTotalTravelCost(0);
        ch.setFitness(Double.POSITIVE_INFINITY);
        double travelCost = 0;
        double[] highestAndTotalTardiness = new double[2];
        int[] routeEndPoint = new int[allCaregivers.length];
        double[] routesCurrentTime = new double[allCaregivers.length];

        Set<Integer> track = trackHolder.get(); // Initial capacity for small patient sets
        List<Integer>[] genes = ch.getGenes();
        for (int i = 0; i < genes.length; i++) {
            List<Integer> route = genes[i];
            for (int j = 0; j < route.size(); j++) {
                if(j == 0){
                    int nextIndex = route.get(j)+1;
                    travelCost += distances[0][nextIndex];
                    if(route.size() == 1){
                        travelCost += distances[nextIndex][0];
                    }
                }else if(j == route.size()-1){
                    int prevIndex = route.get(j-1)+1;
                    int nextIndex = route.get(j)+1;
                    travelCost += distances[prevIndex][nextIndex];
                    travelCost += distances[nextIndex][0];
                }else {
                    int prevIndex = route.get(j-1)+1;
                    int nextIndex = route.get(j)+1;
                    travelCost += distances[prevIndex][nextIndex];
                }
            }
        }

        ch.setTotalTravelCost(travelCost);
        travelCost = 1/3d * travelCost;
        try {
            double solutionCost =0;
            for (int i = 0; i < routeEndPoint.length; i++) {
                List<Integer> route = genes[i];
                int routeEnd = routeEndPoint[i];
                if (routeEnd != -1) {
                    for (int j = routeEnd; j < route.size(); j++) {
                        if (j == 0) {
                            solutionCost = patientIsAssigned(genes, i, -1, route.get(j), travelCost, routesCurrentTime, highestAndTotalTardiness, routeEndPoint, track);
                        } else {
                            solutionCost = patientIsAssigned(genes, i,  route.get(j - 1), route.get(j), travelCost, routesCurrentTime, highestAndTotalTardiness, routeEndPoint, track);
                        }
                        if (solutionCost == Double.POSITIVE_INFINITY) {
                            ch.setFitness(Double.POSITIVE_INFINITY);
                            return;
                        }

                        track.clear();
                    }
                }
            }
            ch.setHighestTardiness(highestAndTotalTardiness[0]);
            ch.setTotalTardiness(highestAndTotalTardiness[1]);
            ch.setFitness(solutionCost);
        }
        finally {
            track.clear();
        }
    }

    private static double patientIsAssigned(List<Integer>[] genes, int route1, int curPatient1, int patient, double totalTravelCost, double[] routesCurrentTime, double[] highestAndTotalTardiness, int[] routeEndPoint, Set<Integer> track) {
        Patient p = allPatients[patient];
        double[] timeWindow = p.getTime_window();
        int currentLocation1 = curPatient1 + 1;
        int nextLocation = patient + 1;

        double arrivalTime1 = routesCurrentTime[route1] + distances[currentLocation1][nextLocation];
        double startTime1 = Math.max(arrivalTime1, timeWindow[0]);

        if (p.getRequired_caregivers().length > 1) {
            if (!track.add(patient)) {
                return Double.POSITIVE_INFINITY;
            }

            int route2 = findOtherCaregiver(patient, route1, genes, totalTravelCost,  routesCurrentTime, highestAndTotalTardiness, routeEndPoint, track);
            if (route2 > allCaregivers.length - 1) {
                return Double.POSITIVE_INFINITY;
            }

            int position2 = routeEndPoint[route2] - 1;
            int curPatient2;
            if (position2 == -1) {
                curPatient2 = -1;
            }else {
                curPatient2 = genes[route2].get(position2);
            }

            if (SwapRoutes(p, route1, route2)) {
                int temp = route1;
                route1 = route2;
                route2 = temp;
                temp = curPatient1;
                curPatient1 = curPatient2;
                curPatient2 = temp;
                currentLocation1 = curPatient1 + 1;
                arrivalTime1 = routesCurrentTime[route1] + distances[currentLocation1][nextLocation];
                startTime1 = Math.max(arrivalTime1, timeWindow[0]);
            }

            int currentLocation2 = curPatient2 + 1;
            double arrivalTime2 = routesCurrentTime[route2] + distances[currentLocation2][nextLocation];
            double startTime2 = Math.max(arrivalTime2, timeWindow[0]);

            if (p.getSynchronization().getType().equals("sequential")) {
                double[] syncDistances = p.getSynchronization().getDistance();
                startTime2 = Math.max(startTime2, startTime1 + syncDistances[0]);
                if (startTime2 - startTime1 > syncDistances[1]) {
                    startTime1 = startTime2 - syncDistances[1];
                }
                double tardiness = Math.max(0, startTime1 - timeWindow[1]);
                highestAndTotalTardiness[0] = Math.max(highestAndTotalTardiness[0], tardiness);
                highestAndTotalTardiness[1] += tardiness;
                tardiness = Math.max(0, startTime2 - timeWindow[1]);
                highestAndTotalTardiness[0] = Math.max(highestAndTotalTardiness[0], tardiness);
                highestAndTotalTardiness[1] += tardiness;
            }else {
                double startTime = Math.max(startTime1, startTime2);
                double tardiness = Math.max(0, startTime - timeWindow[1]);
                highestAndTotalTardiness[0] = Math.max(highestAndTotalTardiness[0], tardiness);
                highestAndTotalTardiness[1] += (2*tardiness);
                startTime1 = startTime2 = startTime;
            }
            routeEndPoint[route1]++;
            routeEndPoint[route2]++;
            routesCurrentTime[route1] = startTime1 + p.getRequired_caregivers()[0].getDuration();
            routesCurrentTime[route2] = startTime2 + p.getRequired_caregivers()[1].getDuration();

        } else {
            double tardiness = Math.max(0, startTime1 - timeWindow[1]);
            highestAndTotalTardiness[0] = Math.max(highestAndTotalTardiness[0], tardiness);
            highestAndTotalTardiness[1] += tardiness;
            routeEndPoint[route1]++;
            routesCurrentTime[route1] = startTime1 + p.getRequired_caregivers()[0].getDuration();
        }
        return totalTravelCost + (1 / 3d * highestAndTotalTardiness[0]) + (1 / 3d * highestAndTotalTardiness[1]);
    }

    private static int findOtherCaregiver(int patient, int route1, List<Integer>[] genes, double totalTravelCost, double[] routesCurrentTime, double[] highestAndTotalTardiness, int[] routeEndPoint, Set<Integer> track) {
        for (int i = 0; i < genes.length; i++) {
            if (i != route1 && genes[i].contains(patient)) {
                List<Integer> route = genes[i];
                int patientPosition = route.indexOf(patient);

                int j = routeEndPoint[i];
                while(routeEndPoint[i] != patientPosition && j < route.size()) {
                    int curPatient;
                    if(j==0)
                        curPatient = -1;
                    else
                        curPatient = route.get(j-1);
                    if(patientIsAssigned(genes,i,curPatient,route.get(j),totalTravelCost,routesCurrentTime,highestAndTotalTardiness,routeEndPoint,track)==Double.POSITIVE_INFINITY)
                        return Integer.MAX_VALUE;
                    j++;
                }
                return i;
            }
        }
        return Integer.MAX_VALUE;
    }
    private static boolean SwapRoutes(Patient p, int route1, int route2) {
        Required_Caregiver[] requiredCaregivers = p.getRequired_caregivers();
        String service1 = requiredCaregivers[0].getService();
        String service2 = requiredCaregivers[1].getService();

        Set<Integer> service1Routes = dataset.getQualifiedCaregiver(service1);
        Set<Integer> service2Routes = dataset.getQualifiedCaregiver(service2);

        return (service1Routes.contains(route2) && !service2Routes.contains(route2)) ||
                (service2Routes.contains(route1) && !service1Routes.contains(route1)) ||
                (service1Routes.contains(route2) && !service1Routes.contains(route1));
    }

    public static boolean patientAssignment(Chromosome ch, int patient, Shift caregiver1, Shift[] routes, int i, Set<Integer> track) {
        Patient p = allPatients[patient];
        double[] timeWindow = p.getTime_window();
        int currentLocation1 = caregiver1.getRoute().isEmpty() ? 0 : caregiver1.getRoute().get(caregiver1.getRoute().size() - 1) + 1;
        int nextLocation = patient + 1;

        double arrivalTime1 = caregiver1.getCurrentTime().get(caregiver1.getCurrentTime().size() - 1) + distances[currentLocation1][nextLocation];
        double startTime1 = Math.max(arrivalTime1, timeWindow[0]);

        if (p.getRequired_caregivers().length > 1) {
            if (!track.add(patient)) { // Combined contains check and add
                return true;
            }

            int index = findSecondCaregiver(patient, i, routes, ch, track);
            if (index > allCaregivers.length - 1) {
                return true;
            }

            Shift caregiver2 = routes[index];
            Required_Caregiver[] requiredCaregivers = p.getRequired_caregivers();

            // Check and potentially swap caregivers based on qualifications
            if (shouldSwapCaregivers(p, caregiver1, caregiver2)) {
                Shift temp = caregiver1;
                caregiver1 = caregiver2;
                caregiver2 = temp;

                // Recalculate after swap
                currentLocation1 = caregiver1.getRoute().isEmpty() ? 0 : caregiver1.getRoute().get(caregiver1.getRoute().size() - 1) + 1;
                arrivalTime1 = caregiver1.getCurrentTime().get(caregiver1.getCurrentTime().size() - 1) + distances[currentLocation1][nextLocation];
                startTime1 = Math.max(arrivalTime1, timeWindow[0]);
            }

            int currentLocation2 = caregiver2.getRoute().isEmpty() ? 0 : caregiver2.getRoute().get(caregiver2.getRoute().size() - 1) + 1;
            double arrivalTime2 = caregiver2.getCurrentTime().get(caregiver2.getCurrentTime().size() - 1) + distances[currentLocation2][nextLocation];
            double startTime2 = Math.max(arrivalTime2, timeWindow[0]);

            processSynchronization(ch, p, caregiver1, caregiver2, requiredCaregivers,
                    startTime1, startTime2, timeWindow[1]);

            double travelCost = distances[currentLocation1][nextLocation] +
                    distances[currentLocation2][nextLocation];
            updateCaregiverRoutes(ch, caregiver1, caregiver2, patient,
                    distances[currentLocation1][nextLocation],
                    distances[currentLocation2][nextLocation],
                    travelCost);
        } else {
            processSingleCaregiver(ch, caregiver1, p, patient, startTime1, timeWindow[1],
                    distances[currentLocation1][nextLocation]);
        }
        return false;
    }

    static void UpdateCost(Chromosome ch) {
        double cost = (1 / 3d * ch.getTotalTravelCost()) + (1 / 3d * ch.getTotalTardiness()) + (1 / 3d * ch.getHighestTardiness());
//        cost = Math.round(cost * 1000.0) / 1000.0;
        ch.setFitness(cost);
    }

    private static int findSecondCaregiver(int pIndex, int route1, Shift[] routes, Chromosome ch, Set<Integer> track) {
        List<Integer>[] genes = ch.getGenes();

        for (int i = 0; i < genes.length; i++) {
            if (i != route1 && genes[i].contains(pIndex)) {
                List<Integer> route = genes[i];
                Shift caregiver = routes[i];
                int patientPosition = route.indexOf(pIndex);

                // Process patients up to the target position
                int j = caregiver.getRoute().size();
                while (caregiver.getRoute().size() != patientPosition && j < route.size()) {
                    int patient = route.get(j);
                    if (patientAssignment(ch, patient, caregiver, routes, i, track))
                        return Integer.MAX_VALUE;
                    j++;
                }
                return i;
            }
        }
        return Integer.MAX_VALUE;
    }

    private static boolean shouldSwapCaregivers(Patient p, Shift caregiver1, Shift caregiver2) {
        Required_Caregiver[] requiredCaregivers = p.getRequired_caregivers();
        String service1 = requiredCaregivers[0].getService();
        String service2 = requiredCaregivers[1].getService();

        Set<Integer> service1Routes = dataset.getQualifiedCaregiver(service1);
        Set<Integer> service2Routes = dataset.getQualifiedCaregiver(service2);

        int caregiver2Id = caregiver2.getCaregiver().getCacheId();
        int caregiver1Id = caregiver1.getCaregiver().getCacheId();

        return (service1Routes.contains(caregiver2Id) && !service2Routes.contains(caregiver2Id)) ||
                (service2Routes.contains(caregiver1Id) && !service1Routes.contains(caregiver1Id)) ||
                (service1Routes.contains(caregiver2Id) && !service1Routes.contains(caregiver1Id));
    }

    private static void processSynchronization(Chromosome ch, Patient p, Shift caregiver1,
                                               Shift caregiver2, Required_Caregiver[] requiredCaregivers, double startTime1,
                                               double startTime2, double timeWindowEnd) {

        double tardiness1, tardiness2;
        if (p.getSynchronization().getType().equals("sequential")) {
            double[] syncDistances = p.getSynchronization().getDistance();
            startTime2 = Math.max(startTime2, startTime1 + syncDistances[0]);
            if (startTime2 - startTime1 > syncDistances[1]) {
                startTime1 = startTime2 - syncDistances[1];
            }
            tardiness1 = Math.max(0, startTime1 - timeWindowEnd);
            tardiness2 = Math.max(0, startTime2 - timeWindowEnd);
        } else {
            double startTime = Math.max(startTime1, startTime2);
            tardiness1 = tardiness2 = Math.max(0, startTime - timeWindowEnd);
            startTime1 = startTime2 = startTime;
        }

        double maxTardiness = Math.max(tardiness1, tardiness2);
        ch.updateTotalTardiness(tardiness1 + tardiness2);
        ch.setHighestTardiness(Math.max(maxTardiness, ch.getHighestTardiness()));

        caregiver1.setCurrentTime(startTime1 + requiredCaregivers[0].getDuration());
        caregiver1.updateTardiness(tardiness1);
        caregiver2.setCurrentTime(startTime2 + requiredCaregivers[1].getDuration());
        caregiver2.updateTardiness(tardiness2);
    }

    private static void updateCaregiverRoutes(Chromosome ch, Shift caregiver1,
                                              Shift caregiver2, int patientId, double cost1, double cost2, double totalCost) {
        ch.updateTotalTravelCost(totalCost);
        caregiver1.updateRoute(patientId);
        caregiver1.updateTravelCost(cost1);
        caregiver2.updateRoute(patientId);
        caregiver2.updateTravelCost(cost2);
    }

    private static void processSingleCaregiver(Chromosome ch, Shift caregiver,
                                               Patient p, int patient, double startTime, double timeWindowEnd, double travelCost) {
        double tardiness = Math.max(0, startTime - timeWindowEnd);
        ch.setHighestTardiness(Math.max(tardiness, ch.getHighestTardiness()));
        ch.updateTotalTardiness(tardiness);
        ch.updateTotalTravelCost(travelCost);

        caregiver.setCurrentTime(startTime + p.getRequired_caregivers()[0].getDuration());
        caregiver.updateRoute(patient);
        caregiver.updateTravelCost(travelCost);
        caregiver.updateTardiness(tardiness);
    }
}
