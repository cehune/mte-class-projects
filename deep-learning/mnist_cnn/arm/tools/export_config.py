"""Converts the best Optuna trial into a config file for the native trainer."""
import optuna

study = optuna.load_study(study_name="mnist_cnn", storage="sqlite:///../results/optuna.db")
p = study.best_params
out = {
    "n_blocks": p["n_blocks"], "base_filters": p["base_filters"],
    "use_batchnorm": int(p["use_batchnorm"]), "gap": int(p["head"] == "gap"),
    "batch_size": p["batch_size"], "onecycle": 1, "dropout": p["dropout"],
    "weight_decay": p["weight_decay"], "label_smoothing": p["label_smoothing"], "lr": p["lr"],
    "rotation": p["rotation"], "translate": p["translate"], "scale_jitter": p["scale_jitter"],
}
with open("configs/best.txt", "w") as f:
    for key, value in out.items():
        f.write(f"{key} {value}\n")
print("trial", study.best_trial.number, "validation accuracy", study.best_value)
print(out)
