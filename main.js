import * as cal from './modules/mycal.js';

// DOM eles
let ratingsFormEles = document.querySelectorAll('.ratings form');
let ratingsContainerEle = document.getElementById('ratings__container');
let ratingsGenBtnEle = document.getElementById('ratings__btn__gen');
let ratingsRemBtnEle = document.getElementById('ratings__btn__rem');
let calPercEle = document.getElementById('cal__perc');
let calBtnEle = document.getElementById('cal__btn');
let anovaEle = document.getElementsByClassName('anova')[0];
let anovaTableTdEles = anovaEle.querySelectorAll('#anova__table [id]');
let fvaluesEle = document.getElementsByClassName('fvalues')[0];
let fvaluesTableTdEles = fvaluesEle.querySelectorAll('#fvalues__table [id]');
let tukeyEle = document.getElementsByClassName('tukey')[0];
let noTukeyEle = document.getElementsByClassName('no-tukey')[0];
let suitableSampleEle = document.getElementsByClassName('suitable-sample')[0];

const ALPHABET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.split('');

let calPerc,
    judges,
    samples,
    sampleTotalSum,
    sampleTotalSquaredSum,
    judgeRatingTotalSquaredSum,
    sampleRatingSquaredSum,
    dfSamples,
    dfError,
    msError,
    fSamples,
    criticalF,
    samplesMean,
    biggestMean;

function genRatingsTable() {
    // existance check
    if (ratingsContainerEle.hasChildNodes()) {
        alert('A table already exists!');
        return;
    }

    let table = document.createElement('table');
    table.setAttribute('id', 'ratings__table');

    judges = Number(document.getElementById('judges__input').value);
    samples = Number(document.getElementById('samples__input').value);

    // generate table
    for (let i = 0; i <= judges; i++) {
        let tr = document.createElement('tr');
        for (let j = 0; j <= samples; j++) {
            if (i === 0) {
                let th = document.createElement('th');
                let thTxt = j === 0 ? 'Judges' : `Sample ${ALPHABET[j - 1]}`;
                th.appendChild(document.createTextNode(thTxt))
                tr.appendChild(th);;
            } else {
                let td = document.createElement('td');
                if (j === 0) {
                    td.appendChild(document.createTextNode(i));
                } else {
                    let input = document.createElement('input');
                    let attributes = {'type':'number', 'min':'0', 'max':'9', 'required':''};
                    for (let key in attributes) {
                        input.setAttribute(key, attributes[key]);
                    }
                    td.appendChild(input);
                }
                
                tr.appendChild(td);
            }
        }
        
        table.appendChild(tr);
    }

    ratingsContainerEle.appendChild(table);

}


function validInput() {
    let invalid = 0;
    ratingsFormEles.forEach((form) => {
        invalid = form.querySelectorAll(':invalid').length;
    });
    return invalid == 0;
}


function addTotal() {
    let ratingsTable = document.getElementById('ratings__table'),
        headingRow = ratingsTable.firstElementChild,
        totalRow = document.createElement('tr'),
        cols = ratingsTable.firstElementChild.childElementCount,
        rows = ratingsTable.childElementCount;

    sampleTotalSum = 
    sampleTotalSquaredSum = 
    judgeRatingTotalSquaredSum = 
    sampleRatingSquaredSum = 0;
    samplesMean = [];

    totalRow.setAttribute('id', 'ratings__table__totalrow');

    // add total column
    {
        let th = document.createElement('th');
        th.appendChild(document.createTextNode('Total'));
        headingRow.appendChild(th);
    }

    // add total row
    {
        let th = document.createElement('th');
        th.appendChild(document.createTextNode('Total'));
        totalRow.appendChild(th);
        ratingsTable.appendChild(totalRow)
    }

    // add horizontal totals
    for (let i = 1; i <= rows - 1; i++) {
        let total = 0,
            currentRow = ratingsTable.childNodes[i],
            th = document.createElement('th');
        
        for (let j = 1; j <= cols - 1; j++) {
            let currentCellInput = currentRow.childNodes[j].firstElementChild;
            total += Number(currentCellInput.value);
        }

        th.appendChild(document.createTextNode(total));
        currentRow.appendChild(th);

        judgeRatingTotalSquaredSum += total ** 2;
    }

    // add vertical totals
    for (let i = 1; i <= cols - 1; i++) {
        let total = 0,
            th = document.createElement('th');

        for (let j = 1; j <= rows - 1; j++) {
            let currentCell = ratingsTable.childNodes[j].childNodes[i];
            let currentCellInput = currentCell.firstElementChild;
            total += Number(currentCellInput.value)
            sampleRatingSquaredSum += Number(currentCellInput.value) ** 2;
        }
       
        th.appendChild(document.createTextNode(total));
        totalRow.appendChild(th);
        sampleTotalSum += total;
        sampleTotalSquaredSum += total ** 2;
        samplesMean.push([`${ALPHABET[i - 1]}`, cal.round(total / judges)]);
        
    }

    // add sampleTotalSum
    {
        let th = document.createElement('th');
        th.setAttribute('id', 'table__sampleTotalSum');
        th.appendChild(document.createTextNode(sampleTotalSum));
        totalRow.appendChild(th);
    }

}


function genAnovaTable() {
    let cf = cal.cf(sampleTotalSum, judges, samples),
        ssSamples = cal.ssSamples(sampleTotalSquaredSum, judges, cf),
        ssJudges = cal.ssJudges(judgeRatingTotalSquaredSum, samples, cf),
        ssTotal = cal.ssTotal(sampleRatingSquaredSum, cf),
        ssError = cal.ssError(ssTotal, ssSamples, ssJudges),
        dfJudges = cal.dfJudges(judges),
        dfTotal = cal.dfTotal(judges, samples);
    
    dfSamples = cal.dfSamples(samples);
    dfError = cal.dfError(dfTotal, dfSamples, dfJudges);
    msError = cal.msError(ssError, dfError);

    let msSamples = cal.msSamples(ssSamples, dfSamples),
        msJudges = cal.msJudges(ssJudges, dfJudges),
        fJudges = cal.fJudges(msJudges, msError);

    fSamples = cal.fSamples(msSamples, msError);

    let values = [
        dfSamples,
        ssSamples,
        msSamples,
        fSamples,
        dfJudges,
        ssJudges,
        msJudges,
        fJudges,
        dfError,
        ssError,
        msError,
        dfTotal,
        ssTotal
    ];

    // add values
    anovaTableTdEles.forEach((td, index) => {
        let val = isNaN(values[index]) ? 'Undefined' : values[index];
        td.innerHTML = val;
    });

    // show table
    if (anovaEle.style.display != 'block') {
        anovaEle.style.display = 'block';
    }

}


function genFvaluesTable() {
    let fvalueTxt;

    // get calPerc
    calPerc = Number(calPercEle.options[calPercEle.selectedIndex].value);

    if (calPerc == 0.01) {
        fvalueTxt = 'F-value at 1% significance level';
        criticalF = cal.getCriticalF(0.01, dfSamples, dfError);
    } else {
        fvalueTxt = 'F-value at 5% significance level';
        criticalF = cal.getCriticalF(0.05, dfSamples, dfError);
    }

    fvaluesTableTdEles[0].innerHTML = isNaN(fSamples) ? 'Undefined' : fSamples;
    fvaluesTableTdEles[1].innerHTML = fvalueTxt;
    fvaluesTableTdEles[2].innerHTML = criticalF;

    // show table
    if (fvaluesEle.style.display != 'block') {
        fvaluesEle.style.display = 'block';
    }

}


function fvaluesDiffer() {
    return fSamples >= criticalF;
}


function genTukey() {
    // hide notukey ele
    if (noTukeyEle.style.display == 'block') {
        noTukeyEle.style.display = 'none';
    }

    // remove existing ol
    if (tukeyEle.lastElementChild.firstElementChild != null) {
        tukeyEle.lastElementChild.removeChild(tukeyEle.lastElementChild.firstElementChild);
    }
    
    let stdErr = cal.stdErr(msError, judges);
    let lsd = cal.lsd(calPerc, samples, dfError, stdErr);
    let comparedMeans = {};
    let ol = document.createElement('ol');

    // show lsd value
    tukeyEle.getElementsByClassName('lsd')[0].innerHTML = lsd;

    // sort sample means in descending order
    samplesMean.sort((a, b) => {
        return b[1] - a[1];
    });


    biggestMean = samplesMean[0];

    // compare sample means
    while (samplesMean.length > 1) {
        let sample = samplesMean.splice(0, 1);
        samplesMean.forEach((e) => {
            let key = `${sample[0][0]} - ${e[0]}`;
            let val = cal.round(sample[0][1] - e[1]).toFixed(2);
            comparedMeans[key] = val;
        });
    }

    for (let k in comparedMeans) {
        let compTxt = comparedMeans[k] > lsd ? `(>${lsd}) ${'\u2705'}` : `(<${lsd}) ${'\u274C'}`;
        let liTxt = `${k} = ${comparedMeans[k]} ${compTxt}`;
        let li = document.createElement('li');
        li.appendChild(document.createTextNode(liTxt));
        ol.appendChild(li);
    }

    tukeyEle.getElementsByClassName('comparison')[0].appendChild(ol);

    // show tukey ele
    if (tukeyEle.style.display !== 'flex') {
        tukeyEle.style.display = 'flex';
    }

}

function suitableSample() {
    let span = suitableSampleEle.getElementsByTagName('span')[0];
    span.innerHTML = biggestMean[0];

    if (suitableSampleEle.style.display != 'flex') {
        suitableSampleEle.style.display = 'flex';
    }
}


function noTukey() {
    // hide tukey ele
    if (tukeyEle.style.display != 'none') {
        tukeyEle.style.display = 'none';
    }

    // hide suitable sample ele
    if (suitableSampleEle.style.display != 'none') {
        suitableSampleEle.style.display = 'none';
    } 

    let valueEles = noTukeyEle.querySelectorAll('[id]');

    valueEles[0].innerHTML = fSamples;
    valueEles[1].innerHTML = criticalF;
    valueEles[2].innerHTML = calPerc == 0.01 ? '1%' : '5%';


    // show notukey ele
    if (noTukeyEle.style.display !== 'block') {
        noTukeyEle.style.display = 'block';
    }
}


function remTables() {
    if (ratingsContainerEle.hasChildNodes()) {
        ratingsContainerEle.removeChild(ratingsContainerEle.firstElementChild);
        anovaEle.style.display = 
        fvaluesEle.style.display = 
        noTukeyEle.style.display = 
        tukeyEle.style.display  = 
        suitableSampleEle.style.display = 'none';
    }
}


function remTotal() {
    let ratingsTable = document.getElementById('ratings__table');
    if (ratingsTable.querySelector('#ratings__table__totalrow') != null) {
        // remove ratings table's total row
        ratingsTable.removeChild(ratingsTable.lastChild);
        // remove ratings table's total col
        ratingsTable.childNodes.forEach((tr) => {
            tr.removeChild(tr.lastChild);
    });
    }
}


ratingsGenBtnEle.addEventListener('click', () => {
    if (validInput()) {
        genRatingsTable();
    }
})

ratingsRemBtnEle.addEventListener('click', remTables)

calBtnEle.addEventListener('click', () => {
    if (!ratingsContainerEle.hasChildNodes()) {
        alert('Please generate a ratings table first!');
        return;
    } else {
        remTotal();
        addTotal();
        genAnovaTable();
        genFvaluesTable();
        if (fvaluesDiffer()) {
            genTukey();
            suitableSample();
        } else {
            noTukey();
        }
    }
});
