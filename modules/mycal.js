import * as jstat from './myjstat.js';

// round num to two decimal places
export function round(num) {
    return Math.round(num * 100) / 100;
}


export function cf(total, judges, samples) {
    let result = (total ** 2) / (judges * samples);
    return result.toFixed(3);
}


export function ssSamples(sampleTotalSquaredSum, judges, cf) {
    let result = (sampleTotalSquaredSum / judges) - cf;
    return result.toFixed(3);
}


export function ssJudges(judgeRatingTotalSquaredSum, samples, cf) {
    let result = (judgeRatingTotalSquaredSum / samples) - cf;
    return result.toFixed(3);
}


export function ssTotal(sampleRatingSquaredSum, cf) {
    let result = sampleRatingSquaredSum - cf;
    return result.toFixed(3);
}


export function ssError(ssTotal, ssSamples) {
    let result = ssTotal - ssSamples;
    return result.toFixed(3);
}


export function dfSamples(samples) {
    let result = samples - 1;
    return result.toFixed(3);
}


export function dfJudges(judges) {
    let result = judges - 1;
    return result.toFixed(3);
}


export function dfTotal(judges, samples) {
    let result = (judges * samples) - 1;
    return result.toFixed(3);
}


export function dfError(dfTotal, dfSamples) {
    let result = dfTotal - dfSamples;
    return result.toFixed(3);
}


export function msSamples(ssSamples, dfSamples) {
    let result = ssSamples / dfSamples;
    return result.toFixed(3);
}


export function msJudges(ssJudges, dfJudges) {
    let result = ssJudges / dfJudges;
    return result.toFixed(3);
}


export function msError(ssError, dfError) {
    let result = ssError / dfError;
    return result.toFixed(3);
}


export function fSamples(msSamples, msError) {
    let result = msSamples / msError;
    return result.toFixed(3);
}


export function fJudges(msJudges, msError) {
    let result = msJudges / msError;
    return result.toFixed(3);
}

export function stdErr(msError, Judges) {
    let result = Math.sqrt(msError / Judges);
    return result.toFixed(3);
}

export function lsd(calPerc, samples, dfError, stdErr) {
    let result = getCriticalQ(calPerc, samples, dfError) * stdErr;
    return result.toFixed(3);
}

export function getCriticalF(sigLvl, deNum, deDeno) {
    return round(jstat.centralFInv(1 - sigLvl, deNum, deDeno));
}


export function getCriticalQ(sigLvl, deNum, deDeno) {
    return round(jstat.tukeyInv(1 - sigLvl, deNum, deDeno));
}

