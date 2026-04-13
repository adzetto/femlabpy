from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from datetime import date
from pathlib import Path
import subprocess

REPO_ROOT = Path(__file__).resolve().parent.parent
GITHUB_REPO = "https://github.com/adzetto/femlabpy"
PYPI_PROJECT = "https://pypi.org/project/femlabpy"


@dataclass(frozen=True)
class CommitEntry:
    sha: str
    commit_date: date
    subject: str
    labels: tuple[str, ...]


@dataclass(frozen=True)
class ReleaseEntry:
    tag: str
    release_date: date
    previous_tag: str | None
    commits: tuple[CommitEntry, ...]

    @property
    def version(self) -> str:
        return self.tag.removeprefix("v")

    @property
    def commit_count(self) -> int:
        return len(self.commits)

    @property
    def compare_range(self) -> str:
        if self.previous_tag is None:
            return f"start..{self.tag}"
        return f"{self.previous_tag}..{self.tag}"

    @property
    def compare_url(self) -> str | None:
        if self.previous_tag is None:
            return None
        return f"{GITHUB_REPO}/compare/{self.previous_tag}...{self.tag}"

    @property
    def release_url(self) -> str:
        return f"{GITHUB_REPO}/releases/tag/{self.tag}"

    @property
    def tag_url(self) -> str:
        return f"{GITHUB_REPO}/tree/{self.tag}"

    @property
    def pypi_url(self) -> str:
        return f"{PYPI_PROJECT}/{self.version}/"

    @property
    def day_span(self) -> int | None:
        if self.previous_tag is None:
            return None
        previous_date = _tag_date(self.previous_tag)
        return max((self.release_date - previous_date).days, 0)

    @property
    def change_density(self) -> float | None:
        if self.day_span is None:
            return None
        return self.commit_count / max(self.day_span, 1)

    @property
    def label_counts(self) -> Counter[str]:
        counts: Counter[str] = Counter()
        for commit in self.commits:
            counts.update(commit.labels)
        return counts


def _git(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _tag_date(tag: str) -> date:
    return date.fromisoformat(_git("log", "-1", "--format=%cs", tag))


def _classify_subject(subject: str) -> tuple[str, ...]:
    text = subject.lower()
    labels: list[str] = []
    rules = {
        "Docs": (
            "docs",
            "docstring",
            "readme",
            "manual",
            "website",
            "tutorial",
            "api",
            "help()",
            "help ",
        ),
        "Fixes": (
            "fix",
            "bug",
            "compat",
            "pager",
            "plot",
            "elastic",
            "gmsh",
            "regression",
        ),
        "Features": (
            "add",
            "support",
            "implement",
            "introduce",
            "new",
        ),
        "Release": (
            "release",
            "pypi",
            "version",
            "tag",
        ),
        "Maintenance": (
            "test",
            "ci",
            "workflow",
            "build",
            "cleanup",
            "clean up",
            "refactor",
        ),
    }
    for label, keywords in rules.items():
        if any(keyword in text for keyword in keywords):
            labels.append(label)
    if not labels:
        labels.append("Maintenance")
    return tuple(labels)


def _commits_for_range(previous_tag: str | None, tag: str) -> tuple[CommitEntry, ...]:
    args = ["log", "--reverse", "--format=%H%x09%cs%x09%s"]
    if previous_tag is None:
        args.append(tag)
    else:
        args.append(f"{previous_tag}..{tag}")
    raw = _git(*args)
    commits: list[CommitEntry] = []
    if not raw:
        return tuple()
    for line in raw.splitlines():
        sha, commit_date, subject = line.split("\t", 2)
        commits.append(
            CommitEntry(
                sha=sha,
                commit_date=date.fromisoformat(commit_date),
                subject=subject,
                labels=_classify_subject(subject),
            )
        )
    return tuple(commits)


def collect_releases() -> tuple[ReleaseEntry, ...]:
    tags = [tag for tag in _git("tag", "--sort=v:refname").splitlines() if tag.startswith("v")]
    releases: list[ReleaseEntry] = []
    previous_tag: str | None = None
    for tag in tags:
        releases.append(
            ReleaseEntry(
                tag=tag,
                release_date=_tag_date(tag),
                previous_tag=previous_tag,
                commits=_commits_for_range(previous_tag, tag),
            )
        )
        previous_tag = tag
    return tuple(releases)


def _format_timeline_table(releases: tuple[ReleaseEntry, ...]) -> list[str]:
    lines = [
        "| Version | Date | Range | Commits | Days Since Previous | Commits / Day | Links |",
        "| --- | --- | --- | ---: | ---: | ---: | --- |",
    ]
    for release in reversed(releases):
        span = "N/A" if release.day_span is None else str(release.day_span)
        density = "N/A" if release.change_density is None else f"{release.change_density:.2f}"
        if release.compare_url is None:
            links = (
                f"[tag]({release.tag_url}) · [release]({release.release_url}) · "
                f"[PyPI]({release.pypi_url})"
            )
        else:
            links = (
                f"[tag]({release.tag_url}) · [release]({release.release_url}) · "
                f"[compare]({release.compare_url}) · [PyPI]({release.pypi_url})"
            )
        lines.append(
            f"| `{release.tag}` | {release.release_date.isoformat()} | "
            f"`{release.compare_range}` | {release.commit_count} | {span} | {density} | {links} |"
        )
    return lines


def _dominant_labels(release: ReleaseEntry) -> tuple[str, ...]:
    ranked = [label for label, _ in release.label_counts.most_common(3)]
    return tuple(ranked)


def _release_characterization(release: ReleaseEntry) -> str:
    labels = _dominant_labels(release)
    if not labels:
        return "No classified commits were available for this tagged range."
    if len(labels) == 1:
        return f"This release is dominated by `{labels[0]}` work."
    if len(labels) == 2:
        return f"This release is primarily a `{labels[0]}` + `{labels[1]}` update."
    return (
        f"This release combines `{labels[0]}`, `{labels[1]}`, and `{labels[2]}` work "
        "in the tagged commit range."
    )


def _representative_subjects(release: ReleaseEntry, limit: int = 3) -> tuple[str, ...]:
    subjects: list[str] = []
    for commit in release.commits:
        text = commit.subject.strip()
        if "release" in text.lower() and len(release.commits) > 1:
            continue
        subjects.append(text)
        if len(subjects) == limit:
            break
    if not subjects:
        subjects = [commit.subject.strip() for commit in release.commits[:limit]]
    return tuple(subjects)


def _format_label_vector(release: ReleaseEntry) -> str:
    parts = [f"{label}={count}" for label, count in release.label_counts.most_common()]
    if not parts:
        return "`{}`"
    return "`{" + ", ".join(parts) + "}`"


def _format_release_section(release: ReleaseEntry) -> list[str]:
    lines = [
        f"## {release.tag}",
        "",
        f"- Release date: `{release.release_date.isoformat()}`",
        f"- Git range: `{release.compare_range}`",
        f"- Tagged commit count: `{release.commit_count}`",
        f"- Distribution links: [GitHub release]({release.release_url}) and [PyPI {release.version}]({release.pypi_url})",
    ]
    if release.change_density is None:
        lines.append("- Change density: first tagged release in the current version history")
    else:
        lines.append(
            "- Change density: "
            f"`{release.change_density:.2f}` commits/day using `max(1, \\Delta_r)`"
        )
    lines.append("")
    lines.append("### Programmatic Interpretation")
    lines.append("")
    lines.append(_release_characterization(release))
    lines.append("")
    lines.append(f"- Label vector: {_format_label_vector(release)}")
    lines.append("- Representative changes inferred from commit subjects:")
    for subject in _representative_subjects(release):
        lines.append(f"  - {subject}")
    lines.append("")
    lines.append("### Programmatic Summary")
    lines.append("")
    for label, count in release.label_counts.most_common():
        lines.append(f"- {label}: `{count}` tagged commit(s)")
    lines.append("")
    lines.append("### Commit Subjects")
    lines.append("")
    lines.append("<details>")
    lines.append("<summary>Expand raw commit list</summary>")
    lines.append("")
    for commit in release.commits:
        short_sha = commit.sha[:7]
        labels = ", ".join(commit.labels)
        lines.append(
            f"- [`{short_sha}`]({GITHUB_REPO}/commit/{commit.sha}) "
            f"{commit.subject}  "
            f"`{commit.commit_date.isoformat()}` · `{labels}`"
        )
    lines.append("")
    lines.append("</details>")
    lines.append("")
    return lines


def build_release_history_markdown() -> str:
    releases = collect_releases()
    if not releases:
        return "\n".join(
            [
                "## Release History Unavailable",
                "",
                "Git tags were not available when the docs were built, so the release",
                "history page could not be generated from version control.",
                "",
            ]
        )

    current_tag = releases[-1].tag
    lines = [
        "<!-- Generated by docs/_release_history.py; do not edit by hand. -->",
        "",
        "## Release Model",
        "",
        f"The current tagged release is `{current_tag}`. This page is generated from the",
        "ordered git tag set and the commit subjects contained in each tagged range.",
        "",
        "For release `r` with tag `\\tau_r`, previous tag `\\tau_{r-1}`, and commit set",
        r"`\mathcal{C}_r`, the page computes:",
        "",
        "$$",
        r"\mathcal{C}_r = \{ c \mid c \in (\tau_{r-1}, \tau_r] \}",
        "$$",
        "",
        "$$",
        r"N_r = |\mathcal{C}_r|, \qquad",
        r"\Delta_r = \max\left(1, d_r - d_{r-1}\right), \qquad",
        r"\rho_r = \frac{N_r}{\Delta_r}",
        "$$",
        "",
        "$$",
        r"L_r(\ell) = \sum_{c \in \mathcal{C}_r} \mathbf{1}[\ell \in \lambda(c)]",
        "$$",
        "",
        "where `N_r` is the tagged commit count, `\\Delta_r` is the elapsed day span",
        "between releases with a one-day floor, and `\\rho_r` is the resulting commit",
        r"density. The label score `L_r(\ell)` counts how often a release is associated",
        r"with a classification label `\ell` after keyword bucketing commit subjects into",
        "Docs, Fixes, Features, Release, and Maintenance groups.",
        "",
        "## Timeline",
        "",
        *_format_timeline_table(releases),
        "",
    ]
    for release in reversed(releases):
        lines.extend(_format_release_section(release))
    return "\n".join(lines) + "\n"


def generate_release_history(output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(build_release_history_markdown(), encoding="utf-8")


if __name__ == "__main__":
    generate_release_history(REPO_ROOT / "docs" / "_generated" / "releases_body.md")
