"""Build quizzes.json and quiz.html (self-contained interactive quiz app)."""
import json
import os

from quiz_data import BOOK, SOFTWARE

HERE = os.path.dirname(os.path.abspath(__file__))


def pack(d, group):
    out = []
    for title, qs in d.items():
        out.append(dict(group=group, title=title, questions=[
            dict(q=q, options=o, answer=a, explanation=e) for q, o, a, e in qs]))
    return out


data = pack(BOOK, "Book: Mattuck, A Guide to Feynman Diagrams") + pack(SOFTWARE, "Software: feynlab")
json.dump(data, open(os.path.join(HERE, "quizzes.json"), "w"), indent=1, ensure_ascii=False)
nq = sum(len(x["questions"]) for x in data)

HTML = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Feynman Diagrams & feynlab — Quizzes</title>
<style>
:root{--bg:#faf9f6;--card:#fff;--ink:#1f1f1e;--muted:#6b6a64;--line:#e6e5e0;--acc:#2a78d6;--ok:#0f7b3e;--okbg:#e6f4ea;--bad:#b42318;--badbg:#fdecea}
@media (prefers-color-scheme:dark){:root{--bg:#1a1a19;--card:#232322;--ink:#f0efec;--muted:#a8a79f;--line:#383836;--acc:#3987e5;--ok:#5cc489;--okbg:#17301f;--bad:#f2867c;--badbg:#3a1d1a}}
*{box-sizing:border-box}body{margin:0;font:15px/1.55 system-ui,-apple-system,Segoe UI,sans-serif;background:var(--bg);color:var(--ink)}
.wrap{display:grid;grid-template-columns:290px 1fr;min-height:100vh}
nav{border-right:1px solid var(--line);padding:20px 14px;overflow:auto;max-height:100vh;position:sticky;top:0}
nav h1{font-size:17px;margin:0 0 4px}nav p{color:var(--muted);font-size:12.5px;margin:0 0 14px}
nav h2{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:16px 6px 6px}
nav button{display:flex;justify-content:space-between;width:100%;text-align:left;background:none;border:0;border-radius:8px;padding:7px 9px;color:var(--ink);font:inherit;font-size:13.5px;cursor:pointer}
nav button:hover{background:var(--line)}nav button.on{background:var(--acc);color:#fff}
nav .sc{font-variant-numeric:tabular-nums;opacity:.75;font-size:12px}
main{padding:28px clamp(16px,4vw,56px);max-width:860px}
main h2{margin:0 0 4px;font-size:22px}.sub{color:var(--muted);margin:0 0 22px}
.q{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:0 0 16px}
.q h3{font-size:15.5px;margin:0 0 12px;font-weight:600}.q h3 span{color:var(--muted);font-weight:500;margin-right:6px}
.opt{display:block;width:100%;text-align:left;background:none;border:1px solid var(--line);border-radius:8px;padding:9px 12px;margin:6px 0;color:var(--ink);font:inherit;cursor:pointer}
.opt:hover:not([disabled]){border-color:var(--acc)}
.opt.ok{background:var(--okbg);border-color:var(--ok)}.opt.bad{background:var(--badbg);border-color:var(--bad)}
.opt[disabled]{cursor:default}
.ex{margin-top:10px;font-size:14px;color:var(--muted)}.ex b{color:var(--ink)}
.bar{display:flex;gap:10px;align-items:center;margin:8px 0 22px}
.bar button{background:var(--acc);color:#fff;border:0;border-radius:8px;padding:8px 14px;font:inherit;cursor:pointer}
.bar .ghost{background:none;color:var(--acc);border:1px solid var(--acc)}
.score{font-variant-numeric:tabular-nums;font-weight:600}
@media (max-width:760px){.wrap{grid-template-columns:1fr}nav{position:static;max-height:none;border-right:0;border-bottom:1px solid var(--line)}}
/* gs-watermark */body::before{content:"\00A9  Gokaran Shukla";position:fixed;left:50%;top:50%;transform:translate(-50%,-50%) rotate(-30deg);font:700 11vw/1 system-ui,sans-serif;color:var(--ink);opacity:.05;white-space:nowrap;pointer-events:none;z-index:0}.wrap{position:relative;z-index:1}.gs-foot{position:fixed;bottom:6px;right:12px;font-size:11px;color:var(--muted);opacity:.85;z-index:2;pointer-events:none}
</style></head><body><div class="gs-foot">© 2026 Gokaran Shukla · github.com/Gokaran-Shukla/feynlab · MIT License</div><div class="wrap"><nav id="nav"></nav><main id="main"></main></div>
<script>
const DATA = __DATA__;
let state = {}; try { state = JSON.parse(localStorage.getItem('feynlabQuiz') || '{}'); } catch(e) {}
function save(){ try { localStorage.setItem('feynlabQuiz', JSON.stringify(state)); } catch(e) {} }
let cur = 0;
function score(i){ const s = state[i] || {}; let c = 0, n = 0; DATA[i].questions.forEach((q,j)=>{ if(s[j]!==undefined){ n++; if(s[j]===q.answer) c++; }}); return [c,n]; }
function nav(){
  const el = document.getElementById('nav'); let html = '<h1>Feynman Diagrams Quizzes</h1><p>__NQ__ questions · book chapters + feynlab tutorials · answers saved in this browser</p>';
  let g = null;
  DATA.forEach((qz,i)=>{ if(qz.group!==g){ g=qz.group; html += `<h2>${g}</h2>`; }
    const [c,n] = score(i); html += `<button class="${i===cur?'on':''}" onclick="show(${i})"><span>${qz.title}</span><span class="sc">${n?c+'/'+qz.questions.length:''}</span></button>`; });
  el.innerHTML = html;
}
function show(i){
  cur = i; const qz = DATA[i]; const s = state[i] || {}; const [c,n] = score(i);
  let html = `<h2>${qz.title}</h2><p class="sub">${qz.group}</p><div class="bar"><span class="score">Score: ${c} / ${qz.questions.length}${n<qz.questions.length?' · '+(qz.questions.length-n)+' unanswered':''}</span><button class="ghost" onclick="reset(${i})">Reset</button>${i<DATA.length-1?`<button onclick="show(${i+1})">Next quiz →</button>`:''}</div>`;
  qz.questions.forEach((q,j)=>{
    const done = s[j]!==undefined;
    html += `<div class="q"><h3><span>${j+1}.</span>${q.q}</h3>`;
    q.options.forEach((o,k)=>{ let cls = ''; if(done){ if(k===q.answer) cls='ok'; else if(k===s[j]) cls='bad'; }
      html += `<button class="opt ${cls}" ${done?'disabled':''} onclick="answer(${i},${j},${k})">${String.fromCharCode(65+k)}. ${o}</button>`; });
    if(done) html += `<div class="ex"><b>${s[j]===q.answer?'Correct.':'Not quite — answer '+String.fromCharCode(65+q.answer)+'.'}</b> ${q.explanation||''}</div>`;
    html += '</div>';
  });
  document.getElementById('main').innerHTML = html; nav(); window.scrollTo(0,0);
}
function answer(i,j,k){ state[i] = state[i] || {}; state[i][j] = k; save(); const y = window.scrollY; show(i); window.scrollTo(0,y); }
function reset(i){ delete state[i]; save(); show(i); }
show(0);
</script></body></html>"""
html = HTML.replace("__DATA__", json.dumps(data, ensure_ascii=False)).replace("__NQ__", str(nq))
open(os.path.join(HERE, "quiz.html"), "w").write(html)

# markdown versions (printable, with answer key)
md = ["# Quizzes — Mattuck chapters and feynlab tutorials", "", "<!-- gs-watermark -->",
      "> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.", ""]
for qz in data:
    md += [f"## {qz['title']}", ""]
    for j, q in enumerate(qz["questions"]):
        md.append(f"{j+1}. {q['q']}")
        md += [f"   - {chr(65+k)}. {o}" for k, o in enumerate(q["options"])]
    md += ["", "<details><summary>Answer key</summary>", ""]
    md += [f"{j+1}. **{chr(65+q['answer'])}** — {q['explanation']}" for j, q in enumerate(qz["questions"])]
    md += ["", "</details>", ""]
open(os.path.join(HERE, "QUIZZES.md"), "w").write("\n".join(md))
print(len(data), "quizzes,", nq, "questions")
